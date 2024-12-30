use crate::math_util::*;
use raylib::prelude::*;
use std::fmt::Debug;

// This is only really a useful abstraction for rendering or calculating collisions #[derive(Debug)]
#[derive(Debug)]
pub struct Polygon {
    // These points are absolute
    pub points: Vec<Vector3>,
}

impl Polygon {
    /**
    Points have to be sorted being passed to the polygon constructor
     */
    pub fn new(points: Vec<Vector3>) -> Self {
        Polygon { points }
    }

    pub fn get_edges(&self) -> Vec<(Vector3, Vector3)> {
        let mut edges: Vec<(Vector3, Vector3)> = Vec::with_capacity(self.points.len());
        for i in 0..self.points.len() {
            edges.push((self.points[i], self.points[(i + 1) % self.points.len()]));
        }
        return edges;
    }

    /**
    Calculates the orthogonal planes necessary for the separating axis theorem
    Some of these planes could be the same as ones already checked, so this isn't quite optimal
     */
    pub fn calculate_orthogonal_planes(&self) -> Vec<Plane> {
        let polygon_as_plane = self.convert_to_plane();
        let mut planes: Vec<Plane> = Vec::with_capacity(self.points.len());
        for edge in self.get_edges() {
            let edge_vec = edge.1 - edge.0;
            let new_normal = polygon_as_plane.n.cross(edge_vec);
            planes.push(Plane::from_point_and_normal(edge.0, new_normal));
        }
        return planes;
    }

    pub fn convert_to_plane(&self) -> Plane {
        // Just make it the first point because why not, we just need *a* point
        let p0 = self.points.first().unwrap();

        // Calculate cross product

        // Every polygon should have three points, so should be able to safely retrieve the second
        // and third points in the vector
        let p1 = self.points.get(1).unwrap();
        let p2 = self.points.get(2).unwrap();

        let v1 = Vector3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);

        let v2 = Vector3::new(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);

        let cross = v1.cross(v2);

        Plane::from_point_and_normal(*p0, cross)
    }

    pub fn collides_with(&self, other: &Polygon) -> bool {
        for i in 0..self.points.len() {
            let line_segment = self.get_edges()[i];
            let third_point = self.get_edges()[(i + 1) % self.points.len()].1; // Borrow another point so we can calculate the normal vector to the line
            let u = line_segment.1 - line_segment.0;
            let v = third_point - line_segment.0;
            let polygon_normal = u.cross(v);
            let normal = u.cross(polygon_normal);
            let projection_line = Line3D::from_point_and_parallel_vec(line_segment.0, normal);
            let line_segment1 = projection_line.project_polygon_onto_line(self);
            let line_segment2 = projection_line.project_polygon_onto_line(other);
            if !line_segments_overlap(line_segment1, line_segment2) {
                return false;
            }
        }

        for i in 0..other.points.len() {
            let line_segment = other.get_edges()[i];
            let third_point = other.get_edges()[(i + 1) % self.points.len()].1;
            let u = line_segment.1 - line_segment.0;
            let v = third_point - line_segment.0;
            let polygon_normal = u.cross(v);
            let normal = u.cross(polygon_normal);
            let projection_line = Line3D::from_point_and_parallel_vec(line_segment.0, normal);
            let line_segment1 = projection_line.project_polygon_onto_line(self);
            let line_segment2 = projection_line.project_polygon_onto_line(other);
            if !line_segments_overlap(line_segment1, line_segment2) {
                return false;
            }
        }

        return true;
    }
}

pub trait MeshShape: Debug {
    fn move_by(&mut self, change: Vector3);
    fn get_vertices(&self) -> Vec<Vector3>;
    fn get_polygons(&self) -> Vec<Polygon>;
    fn get_center(&self) -> Vector3;
    fn get_bounding_circle_radius(&self) -> f32;
    fn render(&self, draw_handle: &mut RaylibMode3D<'_, RaylibDrawHandle<'_>>);
}

#[derive(Debug)]
pub struct RectangularPrism {
    // The root is an absolute point and is the bottom left point of the front face
    pub root: Vector3,

    pub length: f32,
    pub width: f32,
    pub height: f32,

    bounding_circle_radius: f32,
}

impl RectangularPrism {
    pub fn new(root: Vector3, length: f32, width: f32, height: f32) -> Self {
        // Calculate the furthest point
        let bounding_circle_radius = Vector3::new(width / 2.0, height / 2.0, length / 2.0).length();
        RectangularPrism {
            root,
            length,
            width,
            height,
            bounding_circle_radius,
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
            Polygon::new(vec![vertices[0], vertices[1], vertices[3], vertices[2]]),
            Polygon::new(vec![vertices[4], vertices[5], vertices[7], vertices[6]]),
            Polygon::new(vec![vertices[0], vertices[1], vertices[5], vertices[4]]),
            Polygon::new(vec![vertices[2], vertices[3], vertices[7], vertices[6]]),
            Polygon::new(vec![vertices[0], vertices[2], vertices[6], vertices[4]]),
            Polygon::new(vec![vertices[1], vertices[3], vertices[7], vertices[5]]),
        ]
    }

    fn get_center(&self) -> Vector3 {
        return Vector3::new(
            self.root.x + self.width / 2.0,
            self.root.y + self.height / 2.0,
            self.root.z + self.length / 2.0,
        );
    }

    fn get_bounding_circle_radius(&self) -> f32 {
        self.bounding_circle_radius
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
            let polygon_points = polygon.points;

            // TODO: For some reason this works, but I think it shouldn't...
            for i in 1..polygon_points.len() - 1 {
                // Draw the front side
                draw_mode.draw_triangle3D(
                    polygon_points[0],
                    polygon_points[i],
                    polygon_points[i + 1],
                    colors[color_index],
                );

                // Draw the back side
                draw_mode.draw_triangle3D(
                    polygon_points[i + 1],
                    polygon_points[i],
                    polygon_points[0],
                    colors[color_index],
                );
            }

            color_index += 1;
            color_index %= colors.len();
        }
    }
}
