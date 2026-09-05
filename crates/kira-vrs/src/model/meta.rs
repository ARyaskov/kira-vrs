//! Decorative metadata shared by all VRS entities (gkm-core `Entity` / `Element` properties)
//! plus the `Expression` and `Extension` helper classes.
//!
//! Metadata is stored behind an `Option<Box<Meta>>` on every entity, so an object without
//! metadata costs a single null pointer. None of these fields participate in digest
//! computation.

use crate::model::identifier::Digest;
use crate::model::primitives::Syntax;

/// gkm-core `Entity` properties plus VRS `Ga4ghIdentifiableObject.digest` and
/// `Variation.expressions`.
///
/// Which fields are meaningful depends on the owning class: `digest` is emitted only for
/// identifiable objects and `expressions` only for variation objects. Fields are public
/// because this is a plain, unvalidated data bag.
#[derive(Clone, Debug, Default, PartialEq)]
#[non_exhaustive]
pub struct Meta {
    /// Logical identifier within a system, conventionally the GA4GH computed identifier.
    pub id: Option<String>,
    /// A primary name.
    pub name: Option<String>,
    /// Free-text description.
    pub description: Option<String>,
    /// Alternative names.
    pub aliases: Vec<String>,
    /// Extensions carrying data outside the standard.
    pub extensions: Vec<Extension>,
    /// A previously computed sha512t24u digest (carried, not trusted: digests are always
    /// recomputed from content by [`Identifiable::digest`](crate::digest::Identifiable::digest)).
    pub digest: Option<Digest>,
    /// Nomenclature expressions (HGVS, SPDI, ...) of a variation.
    pub expressions: Vec<Expression>,
}

impl Meta {
    /// An empty metadata block.
    pub fn new() -> Self {
        Self::default()
    }

    /// `true` if every field is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_none()
            && self.name.is_none()
            && self.description.is_none()
            && self.aliases.is_empty()
            && self.extensions.is_empty()
            && self.digest.is_none()
            && self.expressions.is_empty()
    }
}

/// Access to the metadata block of any VRS entity.
pub trait Entity {
    /// The metadata block, if any.
    fn meta(&self) -> Option<&Meta>;

    /// Mutable access to the metadata block, creating it if absent.
    fn meta_mut(&mut self) -> &mut Meta;

    /// The `id` property.
    fn id(&self) -> Option<&str> {
        self.meta().and_then(|m| m.id.as_deref())
    }

    /// The `name` property.
    fn name(&self) -> Option<&str> {
        self.meta().and_then(|m| m.name.as_deref())
    }

    /// The `description` property.
    fn description(&self) -> Option<&str> {
        self.meta().and_then(|m| m.description.as_deref())
    }

    /// The `aliases` property.
    fn aliases(&self) -> &[String] {
        self.meta().map_or(&[], |m| m.aliases.as_slice())
    }

    /// The `extensions` property.
    fn extensions(&self) -> &[Extension] {
        self.meta().map_or(&[], |m| m.extensions.as_slice())
    }

    /// Set the `id` property.
    fn set_id(&mut self, id: impl Into<String>) {
        self.meta_mut().id = Some(id.into());
    }
}

/// Implements [`Entity`] and the `with_*` metadata builders for a struct with a
/// `meta: Option<Box<Meta>>` field.
macro_rules! impl_entity {
    ($ty:ty) => {
        impl $crate::model::Entity for $ty {
            #[inline]
            fn meta(&self) -> Option<&$crate::model::Meta> {
                self.meta.as_deref()
            }
            #[inline]
            fn meta_mut(&mut self) -> &mut $crate::model::Meta {
                self.meta.get_or_insert_with(Default::default)
            }
        }

        impl $ty {
            /// Attach a metadata block (replacing any existing one).
            #[must_use]
            pub fn with_meta(mut self, meta: $crate::model::Meta) -> Self {
                self.meta = if meta.is_empty() {
                    None
                } else {
                    Some(Box::new(meta))
                };
                self
            }

            /// Set the `id` property.
            #[must_use]
            pub fn with_id(mut self, id: impl Into<String>) -> Self {
                $crate::model::Entity::meta_mut(&mut self).id = Some(id.into());
                self
            }

            /// Set the `name` property.
            #[must_use]
            pub fn with_name(mut self, name: impl Into<String>) -> Self {
                $crate::model::Entity::meta_mut(&mut self).name = Some(name.into());
                self
            }

            /// Set the `description` property.
            #[must_use]
            pub fn with_description(mut self, description: impl Into<String>) -> Self {
                $crate::model::Entity::meta_mut(&mut self).description = Some(description.into());
                self
            }

            /// Add an alias.
            #[must_use]
            pub fn with_alias(mut self, alias: impl Into<String>) -> Self {
                $crate::model::Entity::meta_mut(&mut self)
                    .aliases
                    .push(alias.into());
                self
            }

            /// Add an extension.
            #[must_use]
            pub fn with_extension(mut self, extension: $crate::model::Extension) -> Self {
                $crate::model::Entity::meta_mut(&mut self)
                    .extensions
                    .push(extension);
                self
            }
        }
    };
}
pub(crate) use impl_entity;

/// Implements the `expressions` accessors for variation classes.
macro_rules! impl_variation_expressions {
    ($ty:ty) => {
        impl $ty {
            /// Nomenclature expressions (HGVS, SPDI, ...) describing this variation.
            pub fn expressions(&self) -> &[$crate::model::Expression] {
                self.meta
                    .as_deref()
                    .map_or(&[], |m| m.expressions.as_slice())
            }

            /// Add a nomenclature expression.
            #[must_use]
            pub fn with_expression(mut self, expression: $crate::model::Expression) -> Self {
                $crate::model::Entity::meta_mut(&mut self)
                    .expressions
                    .push(expression);
                self
            }
        }
    };
}
pub(crate) use impl_variation_expressions;

/// A gkm-core `Extension`: a named value outside the standard model.
///
/// The value is arbitrary JSON by definition, so it is held as a [`serde_json::Value`]; this is
/// the one place in the model where a dynamic JSON value is appropriate.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub struct Extension {
    /// Name indicative of the meaning of the value.
    pub name: String,
    /// The value (any JSON).
    pub value: serde_json::Value,
    /// Description of the meaning or utility of the extension.
    pub description: Option<String>,
    /// Logical identifier.
    pub id: Option<String>,
    /// Nested extensions.
    pub extensions: Vec<Extension>,
}

impl Extension {
    /// Create an extension.
    pub fn new(name: impl Into<String>, value: impl Into<serde_json::Value>) -> Self {
        Self {
            name: name.into(),
            value: value.into(),
            description: None,
            id: None,
            extensions: Vec::new(),
        }
    }

    /// Set the description.
    #[must_use]
    pub fn with_description(mut self, description: impl Into<String>) -> Self {
        self.description = Some(description.into());
        self
    }
}

/// A VRS `Expression`: the variation written in another nomenclature (HGVS, SPDI, ...).
#[derive(Clone, Debug, PartialEq)]
pub struct Expression {
    syntax: Syntax,
    value: String,
    syntax_version: Option<String>,
    id: Option<String>,
    extensions: Vec<Extension>,
}

impl Expression {
    /// Create an expression in the given syntax.
    pub fn new(syntax: Syntax, value: impl Into<String>) -> Self {
        Self {
            syntax,
            value: value.into(),
            syntax_version: None,
            id: None,
            extensions: Vec::new(),
        }
    }

    /// Set the syntax version (important for HGVS, whose syntax has evolved).
    #[must_use]
    pub fn with_syntax_version(mut self, version: impl Into<String>) -> Self {
        self.syntax_version = Some(version.into());
        self
    }

    /// Set the logical identifier.
    #[must_use]
    pub fn with_id(mut self, id: impl Into<String>) -> Self {
        self.id = Some(id.into());
        self
    }

    /// Add an extension.
    #[must_use]
    pub fn with_extension(mut self, extension: Extension) -> Self {
        self.extensions.push(extension);
        self
    }

    /// The nomenclature.
    #[inline]
    pub fn syntax(&self) -> Syntax {
        self.syntax
    }

    /// The expression text.
    #[inline]
    pub fn value(&self) -> &str {
        &self.value
    }

    /// The syntax version, if given.
    #[inline]
    pub fn syntax_version(&self) -> Option<&str> {
        self.syntax_version.as_deref()
    }

    /// The logical identifier, if given.
    #[inline]
    pub fn id(&self) -> Option<&str> {
        self.id.as_deref()
    }

    /// Extensions.
    #[inline]
    pub fn extensions(&self) -> &[Extension] {
        &self.extensions
    }

    /// Construct from all parts (used by deserialization).
    pub(crate) fn from_parts(
        syntax: Syntax,
        value: String,
        syntax_version: Option<String>,
        id: Option<String>,
        extensions: Vec<Extension>,
    ) -> Self {
        Self {
            syntax,
            value,
            syntax_version,
            id,
            extensions,
        }
    }
}
