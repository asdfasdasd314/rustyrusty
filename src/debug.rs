use raylib::prelude::*;
use crate::render::*;

pub fn draw_wireframe(d: &mut RaylibMode3D<'_, RaylibDrawHandle<'_>>, polygons: Vec<Polygon>) {
    for polygon in polygons {
        let points = polygon.points;
        for i in 0..points.len() {
            d.draw_line_3D(points[0], points[i], Color::BLACK);
        }
    }
}
