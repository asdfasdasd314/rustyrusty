use raylib::prelude::{Vector2, Vector3};

/*
To prevent issues with floating point arithmetic, I chose to round numbers in places where the precision matters
Ironically, removing precision will give us precision in calculations because then 5.00000001 == 5.0

This file creates a standardized system of rounding for the entire project

Everything is rounded to the nearest tenth because I believe that is precise enough for this project
*/

// `PRECISION` is the number of decimals of precision
pub const PRECISION: i32 = 1;

pub fn f32_round(x: f32) -> f32 {
    (x * (10.0_f32).powi(PRECISION)).round() / (10.0_f32).powi(PRECISION)
}

pub fn vector3_round(v: Vector3) -> Vector3 {
    Vector3::new(
        f32_round(v.x),
        f32_round(v.y),
        f32_round(v.z),
    )
}

pub fn vector2_round(v: Vector2) -> Vector2 {
    Vector2::new(
        f32_round(v.x),
        f32_round(v.y),
    )
}
