//! Streaming deserialization of `type`-discriminated objects.
//!
//! serde's built-in internally-tagged enums buffer the entire object before dispatching. VRS
//! objects nest several levels deep (`Allele` → `SequenceLocation` → `SequenceReference`), so
//! that would buffer every polymorphic property. Instead, [`TaggedVisitor`] reads map entries
//! until it sees `type`, buffers only the entries that preceded it, and then hands a
//! [`ReplayMap`] — the buffered entries, the `type` entry, then the untouched remainder of the
//! input — to the concrete class's own deserializer.

use std::fmt;

use serde::de::value::{BorrowedStrDeserializer, MapAccessDeserializer, StrDeserializer};
use serde::de::{self, DeserializeSeed, Deserializer, MapAccess, Visitor};

/// A map key or tag string, borrowed from the input when the format allows it.
#[derive(Clone)]
pub(crate) enum KeyStr<'de> {
    Borrowed(&'de str),
    Owned(String),
}

impl KeyStr<'_> {
    fn as_str(&self) -> &str {
        match self {
            Self::Borrowed(s) => s,
            Self::Owned(s) => s,
        }
    }
}

impl<'de> de::Deserialize<'de> for KeyStr<'de> {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = KeyStr<'de>;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str("a string")
            }
            fn visit_borrowed_str<E: de::Error>(self, v: &'de str) -> Result<Self::Value, E> {
                Ok(KeyStr::Borrowed(v))
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Self::Value, E> {
                Ok(KeyStr::Owned(v.to_owned()))
            }
            fn visit_string<E: de::Error>(self, v: String) -> Result<Self::Value, E> {
                Ok(KeyStr::Owned(v))
            }
        }
        deserializer.deserialize_str(V)
    }
}

/// Per-union dispatch on the `type` tag.
pub(crate) trait TagDispatch<'de>: Sized {
    type Output;

    /// The list of accepted `type` values, for error messages.
    const VARIANTS: &'static [&'static str];

    /// Deserialize the object (whose `type` is `tag`) from `deserializer`, which yields the
    /// complete object including its `type` entry.
    fn dispatch<D: Deserializer<'de>>(tag: &str, deserializer: D)
    -> Result<Self::Output, D::Error>;
}

/// Visitor implementing the streaming dispatch. Also accepts strings when `ACCEPT_IRI` is
/// set, for `IriOr`-like unions whose members include a bare IRI.
pub(crate) struct TaggedVisitor<P> {
    pub(crate) expecting: &'static str,
    pub(crate) _dispatch: std::marker::PhantomData<P>,
}

impl<P> TaggedVisitor<P> {
    pub(crate) const fn new(expecting: &'static str) -> Self {
        Self {
            expecting,
            _dispatch: std::marker::PhantomData,
        }
    }
}

impl<'de, P: TagDispatch<'de>> Visitor<'de> for TaggedVisitor<P> {
    type Value = P::Output;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} object with a `type` property", self.expecting)
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
        let mut buffered: Vec<(KeyStr<'de>, serde_json::Value)> = Vec::new();
        loop {
            let Some(key) = map.next_key::<KeyStr<'de>>()? else {
                return Err(de::Error::missing_field("type"));
            };
            if key.as_str() == "type" {
                let tag: KeyStr<'de> = map.next_value()?;
                let tag_text = tag.clone();
                if !P::VARIANTS.contains(&tag_text.as_str()) {
                    return Err(de::Error::unknown_variant(tag_text.as_str(), P::VARIANTS));
                }
                let replay = ReplayMap {
                    buffered: buffered.into_iter(),
                    pending: None,
                    tag: Some(tag),
                    rest: map,
                };
                return P::dispatch(tag_text.as_str(), MapAccessDeserializer::new(replay));
            }
            let value: serde_json::Value = map.next_value()?;
            buffered.push((key, value));
        }
    }
}

enum Pending<'de> {
    Value(serde_json::Value),
    Tag(KeyStr<'de>),
}

/// A `MapAccess` that replays buffered entries and the `type` entry before continuing with
/// the underlying input.
pub(crate) struct ReplayMap<'de, A> {
    buffered: std::vec::IntoIter<(KeyStr<'de>, serde_json::Value)>,
    pending: Option<Pending<'de>>,
    tag: Option<KeyStr<'de>>,
    rest: A,
}

impl<'de, A: MapAccess<'de>> MapAccess<'de> for ReplayMap<'de, A> {
    type Error = A::Error;

    fn next_key_seed<K: DeserializeSeed<'de>>(
        &mut self,
        seed: K,
    ) -> Result<Option<K::Value>, Self::Error> {
        if let Some((key, value)) = self.buffered.next() {
            self.pending = Some(Pending::Value(value));
            return deserialize_key(seed, &key).map(Some);
        }
        if let Some(tag) = self.tag.take() {
            self.pending = Some(Pending::Tag(tag));
            return seed
                .deserialize(BorrowedStrDeserializer::new("type"))
                .map(Some);
        }
        self.rest.next_key_seed(seed)
    }

    fn next_value_seed<V: DeserializeSeed<'de>>(
        &mut self,
        seed: V,
    ) -> Result<V::Value, Self::Error> {
        match self.pending.take() {
            Some(Pending::Value(value)) => seed.deserialize(value).map_err(de::Error::custom),
            Some(Pending::Tag(tag)) => match tag {
                KeyStr::Borrowed(s) => seed.deserialize(BorrowedStrDeserializer::new(s)),
                KeyStr::Owned(s) => seed.deserialize(StrDeserializer::new(&s)),
            },
            None => self.rest.next_value_seed(seed),
        }
    }

    fn size_hint(&self) -> Option<usize> {
        None
    }
}

fn deserialize_key<'de, K: DeserializeSeed<'de>, E: de::Error>(
    seed: K,
    key: &KeyStr<'de>,
) -> Result<K::Value, E> {
    match key {
        KeyStr::Borrowed(s) => seed.deserialize(BorrowedStrDeserializer::new(s)),
        KeyStr::Owned(s) => seed.deserialize(StrDeserializer::new(s)),
    }
}

/// Deserialize a `type`-tagged union from `deserializer`.
pub(crate) fn deserialize_tagged<'de, D, P>(
    deserializer: D,
    expecting: &'static str,
) -> Result<P::Output, D::Error>
where
    D: Deserializer<'de>,
    P: TagDispatch<'de>,
{
    deserializer.deserialize_map(TaggedVisitor::<P>::new(expecting))
}
