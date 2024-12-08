use crate::math_util::*;
use raylib::prelude::*;

pub struct Polygon {
    pub points: Vec<Vector3>,
}

impl Polygon {
    pub fn new(points: Vec<Vector3>) -> Self {
        Polygon { points }
    }

    pub fn as_plane(&self) -> Plane {
        // Just make it the first point because why not, we just need *a* point
        let p0 = self.points.first().unwrap();

        // Calculate cross product

        // Every polygon should have three points, so should be able to safely retrieve the second
        // and third points in the vector
        let p1 = self.points.get(1).unwrap();
        let p2 = self.points.get(2).unwrap();

        let v1 = Vector3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);

        let v2 = Vector3::new(p2.x - p0.x, p2.y - p0.y, p2.z - p0.z);

        let cross = cross_product(v1, v2);

        Plane::from_point_and_normal(*p0, cross)
    }
}

pub trait MeshShape {
    fn get_vertices(&self) -> Vec<Vector3>;
    fn get_polygons(&self) -> Vec<Polygon>;
    fn render(&self, draw_handle: RaylibDrawHandle);
}

pub struct RectangularPrism {
    // The root is the bottom left point of the front face
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

    fn render(&self, draw_handle: RaylibDrawHandle) {
        for polygon in self.get_polygons() {
            // TODO for tomorrow: get these shapes rendering, I want to be able to render by polygons because there is a ton of room to optimize if this path is taken
        }
    }
}

pub struct SolidBody {
    // Both of these positions are measured from the root of the object (whatever that means)
    // The previous position is necessary for determining collisions
    pub position: Vector3,
    pub velocity: Vector3,

    pub mesh: Box<dyn MeshShape>,
}

impl SolidBody {
    pub fn new(mesh: Box<dyn MeshShape>) -> Self {
        SolidBody {
            position: Vector3::new(0.0, 0.0, 0.0),
            velocity: Vector3::new(0.0, 0.0, 0.0),
            mesh,
        }
    }
}

fn check_polygons_collide(polygon1: &Polygon, polygon2: &Polygon) -> bool {
    // This is why calc 3 was useful :)

    let plane1 = polygon1.as_plane();
    let plane2 = polygon2.as_plane();

    // If this is false, then the polygons are parallel and won't collide
    if !check_planes_intersect(&plane1, &plane2) {
        return false;
    }

    // Now we can calculate the line in 3 dimensions where these planes intersect
    let plane_intersection = calculate_line_intersection_between_planes(&plane1, &plane2);

    match plane_intersection {
        // If the polygons are to intersect, then this intersection line should be contained by both polygons so just use polygon 1
        PlaneIntersection::Line(line) => {
            // Use the line's parallel vector to get the bounding angles
            let bounding_angles = calculate_bounding_angles(&line.v, polygon1);

            return check_angles_surround_relative_origin(bounding_angles);
        }
        PlaneIntersection::VerticalLine(_) => {
            // In this case we already have a relative origin, so just create a vector that points straight up and pass to the calculate relative angles function

            // Create the up vector
            let up = Vector3::new(0.0, 1.0, 0.0);

            let bounding_angles = calculate_bounding_angles(&up, polygon1);
            return check_angles_surround_relative_origin(bounding_angles);
        }
        PlaneIntersection::Infinite => {
            // Immediately return true
            return true;
        }
    }
}

// This function checks that given some angles, they surround the relative origin they were calculated from
fn check_angles_surround_relative_origin(angles: Vec<SphericalAngle>) -> bool {
    // TODO For now, brute force
    for i in 0..(angles.len() - 1) {
        for j in (i + 1)..angles.len() {
            // We want to check if the sign flips across the pair in both the theta and phi
            let angle1 = &angles[i];
            let angle2 = &angles[j];

            // Sign should be flipped, so check for that (the product should be less than 0)
            if angle1.theta * angle2.theta <= 0.0 && angle1.phi * angle2.phi <= 0.0 {
                return true;
            }
        }
    }

    return false;
}

fn calculate_bounding_angles(base_vec: &Vector3, polygon: &Polygon) -> Vec<SphericalAngle> {
    let mut bounding_vecs: Vec<Vector3> = Vec::with_capacity(polygon.points.len());
    for vertex in polygon.points.iter() {
        let bounding_vec = Vector3::new(
            vertex.x - base_vec.x,
            vertex.y - base_vec.y,
            vertex.z - base_vec.z,
        );

        bounding_vecs.push(bounding_vec);
    }

    // Find the difference in angles
    let bounding_angles: Vec<SphericalAngle> = bounding_vecs
        .iter()
        .map(|vec| calculate_difference_in_angle(base_vec, vec))
        .collect();

    return bounding_angles;
}
