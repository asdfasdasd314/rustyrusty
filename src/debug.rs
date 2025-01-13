use raylib::prelude::*;
use crate::render::*;

pub fn draw_wireframe(d: &mut RaylibMode3D<'_, RaylibDrawHandle<'_>>, polygons: Vec<Polygon>) {
    for polygon in polygons {
        let points = polygon.points;
        for i in 0..points.len() {
            d.draw_line_3D(Vector3::from(points[0]), Vector3::from(points[i]), Color::BLACK);
        }
    }
}
