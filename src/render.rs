use crate::math_util::*;
use raylib::prelude::*;
use std::fmt::Debug;

// This is only really a useful abstraction for rendering or calculating collisions
#[derive(Debug)]
pub struct Polygon {
    // These points are absolute
    pub points: Vec<Vector3>,
}

impl Polygon {
    pub fn new(mut points: Vec<Vector3>) -> Self {
        // Ensure the points are stored in a cyclical order (sorted)
        sort_points_by_angle_from_centroid(&mut points);
        Polygon { points }
    }

    pub fn get_edges(&self) -> Vec<(Vector3, Vector3)> {
        let mut edges: Vec<(Vector3, Vector3)> = Vec::with_capacity(self.points.len());
        for i in 0..self.points.len() {
            edges.push((self.points[i], self.points[(i + 1) % self.points.len()]));
        }
        return edges;
    }
}

pub fn convert_polygon_to_plane(polygon: &Polygon) -> Plane {
    // Just make it the first point because why not, we just need *a* point
    let p0 = polygon.points.first().unwrap();

    // Calculate cross product

    // Every polygon should have three points, so should be able to safely retrieve the second
    // and third points in the vector
    let p1 = polygon.points.get(1).unwrap();
    let p2 = polygon.points.get(2).unwrap();

    let v1 = Vector3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);

    let v2 = Vector3::new(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);

    let cross = cross_product(&v1, &v2);

    Plane::from_point_and_normal(*p0, cross)
}

pub trait MeshShape: Debug {
    fn move_by(&mut self, change: Vector3);
    fn get_vertices(&self) -> Vec<Vector3>;
    fn get_polygons(&self) -> Vec<Polygon>;
    fn render(&self, draw_handle: &mut RaylibMode3D<'_, RaylibDrawHandle<'_>>);
}

#[derive(Debug)]
pub struct RectangularPrism {
    // The root is an absolute point and is the bottom left point of the front face
    pub root: Vector3,

    pub length: f32,
    pub width: f32,
    pub height: f32,
}

impl RectangularPrism {
    pub fn new(root: Vector3, length: f32, width: f32, height: f32) -> Self {
        RectangularPrism {
            root,
            length,
            width,
            height,
        }
    }
}

impl MeshShape for RectangularPrism {
    fn move_by(&mut self, change: Vector3) {
        self.root += change;
    }

    fn get_vertices(&self) -> Vec<Vector3> {
        vec![
            Vector3::new(self.root.x, self.root.y, self.root.z),
            Vector3::new(self.root.x + self.width, self.root.y, self.root.z),
            Vector3::new(self.root.x, self.root.y + self.height, self.root.z),
            Vector3::new(
                self.root.x + self.width,
                self.root.y + self.height,
                self.root.z,
            ),
            Vector3::new(self.root.x, self.root.y, self.root.z + self.length),
            Vector3::new(
                self.root.x + self.width,
                self.root.y,
                self.root.z + self.length,
            ),
            Vector3::new(
                self.root.x,
                self.root.y + self.height,
                self.root.z + self.length,
            ),
            Vector3::new(
                self.root.x + self.width,
                self.root.y + self.height,
                self.root.z + self.length,
            ),
        ]
    }

    fn get_polygons(&self) -> Vec<Polygon> {
        let vertices = self.get_vertices();

        // Define all of the polygons
        vec![
            Polygon::new(vec![vertices[0], vertices[1], vertices[2], vertices[3]]),
            Polygon::new(vec![vertices[4], vertices[5], vertices[6], vertices[7]]),
            Polygon::new(vec![vertices[0], vertices[1], vertices[4], vertices[5]]),
            Polygon::new(vec![vertices[2], vertices[3], vertices[6], vertices[7]]),
            Polygon::new(vec![vertices[0], vertices[2], vertices[4], vertices[6]]),
            Polygon::new(vec![vertices[1], vertices[3], vertices[5], vertices[7]]),
        ]
    }

    fn render(&self, draw_mode: &mut RaylibMode3D<'_, RaylibDrawHandle<'_>>) {
        let polygons = self.get_polygons();

        let mut color_index = 0;
        let colors = vec![
            Color::RED,
            Color::BLACK,
            Color::YELLOW,
            Color::GREEN,
            Color::BLUE,
            Color::PURPLE,
        ];
        for polygon in polygons {
            let mut polygon_points = polygon.points;
            sort_points_by_angle_from_centroid(&mut polygon_points);

            for i in 0..polygon_points.len() - 2 {
                // Draw the front side
                draw_mode.draw_triangle3D(
                    polygon_points[i],
                    polygon_points[i + 1],
                    polygon_points[i + 2],
                    colors[color_index],
                );

                // Draw the back side
                draw_mode.draw_triangle3D(
                    polygon_points[i + 2],
                    polygon_points[i + 1],
                    polygon_points[i],
                    colors[color_index],
                );
            }

            color_index += 1;
            color_index %= colors.len();
        }
    }
}
