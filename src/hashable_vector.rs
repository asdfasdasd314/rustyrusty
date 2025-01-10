use std::hash::Hash;
use ordered_float::OrderedFloat;
use raylib::prelude::{Vector3, Vector2};

#[derive(Eq, PartialEq, Hash)]
pub struct HashableVector3 {
    pub x: OrderedFloat<f32>,
    pub y: OrderedFloat<f32>,
    pub z: OrderedFloat<f32>,
}

#[derive(Eq, PartialEq, Hash)]
pub struct HashableVector2 {
    pub x: OrderedFloat<f32>,
    pub y: OrderedFloat<f32>,
}

impl HashableVector3 {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self {
            x: OrderedFloat(x),
            y: OrderedFloat(y),
            z: OrderedFloat(z),
        }
    }

    pub fn from_vector3(v: Vector3) -> Self {
        Self {
            x: OrderedFloat(v.x),
            y: OrderedFloat(v.y),
            z: OrderedFloat(v.z),
        }
    }
}

impl HashableVector2 {
    pub fn new(x: f32, y: f32) -> Self {
        Self {
            x: OrderedFloat(x),
            y: OrderedFloat(y),
        }
    }

    pub fn from_vector2(v: Vector2) -> Self {
        Self {
            x: OrderedFloat(v.x),
            y: OrderedFloat(v.y),
        }
    }
}
