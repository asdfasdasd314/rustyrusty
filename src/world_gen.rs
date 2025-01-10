use crate::render::*;
use crate::debug::*;
use crate::math_util::*;
use raylib::prelude::*;
use rand::thread_rng;
use rand::seq::SliceRandom;
use noise::{NoiseFn, Simplex};

#[derive(Debug)]
pub struct GroundMesh {
    pub bottom_left_pos: Vector2,

    // These describe the distance between two points in the heightmap
    pub dx: f32,
    pub dz: f32,

    // Heights of the points on the mesh
    pub h1: f32, // 1st quadrant, upper left
    pub h2: f32, // 2nd quadrant, upper right
    pub h3: f32, // 3rd quadrant, lower right
    pub h4: f32, // 4th quadrant, lower left

    color: raylib::color::Color,
}

impl GroundMesh {
    pub fn new(bottom_left_pos: Vector2, dx: f32, dz: f32, h1: f32, h2: f32, h3: f32, h4: f32) -> Self {
        let colors = vec![
            Color::GREEN,
            Color::BROWN,
            Color::BLACK,
            Color::ORANGE,
            Color::YELLOW,
            Color::PURPLE,
            Color::BLUE,
        ];

        let mut rng = thread_rng();

        let random_color = colors.choose(&mut rng).unwrap();
 
        Self {
            bottom_left_pos,
            dx, 
            dz,
            h1,
            h2,
            h3,
            h4,
            color: *random_color,
        }
    }
}

impl MeshShape for GroundMesh {
    fn move_by(&mut self, change: Vector3) {
        self.bottom_left_pos += Vector2::new(change.x, change.z);
    }

    fn get_vertices(&self) -> Vec<Vector3> {
        return vec![
            Vector3::new(self.bottom_left_pos.x, self.h1, self.bottom_left_pos.y),
            Vector3::new(self.bottom_left_pos.x + self.dx, self.h2, self.bottom_left_pos.y),
            Vector3::new(self.bottom_left_pos.x + self.dx, self.h3, self.bottom_left_pos.y + self.dz),
            Vector3::new(self.bottom_left_pos.x, self.h4, self.bottom_left_pos.y + self.dz),
        ];
    }
    
    fn get_polygons(&self) -> Vec<Polygon> {
        let vertices = self.get_vertices();
        vec![Polygon::new(vertices)]
    }

    fn get_center(&self) -> Vector3 {
        let average_height = (self.h1 + self.h2 + self.h3 + self.h4) / 4.0;
        return Vector3::new(
            self.bottom_left_pos.x + self.dx / 2.0,
            average_height,
            self.bottom_left_pos.y + self.dz / 2.0,
        );
    }

    fn get_bounding_circle_radius(&self) -> f32 {
        let average_height = (self.h1 + self.h2 + self.h3 + self.h4) / 4.0;
        let max_height = f32_max(&[self.h1, self.h2, self.h3, self.h4]);
        (self.dx * self.dx + self.dz * self.dz + max_height - average_height).sqrt()
    }

    fn render(&self, draw_handle: &mut RaylibMode3D<'_, RaylibDrawHandle<'_>>, in_debug: bool) {
        let vertices = self.get_vertices();
       
        for i in 1..vertices.len() - 1 {
            draw_handle.draw_triangle3D(vertices[0], vertices[i], vertices[i + 1], self.color);
            draw_handle.draw_triangle3D(vertices[0], vertices[i + 1], vertices[i], self.color);
        }

        if in_debug {
            draw_wireframe(draw_handle, self.get_polygons());
        }
    }
}

pub fn generate_height_map() -> Vec<Vec<f32>> {
    // Initialize Simplex noise generator
    let simplex = Simplex::new(0);

    // Heightmap dimensions
    let width = 100;
    let height = 100;

    // Generate heightmap using Simplex Noise
    let mut heights = vec![vec![0.0; height]; width];
    for x in 0..width {
        for y in 0..height {
            let noise_value = simplex.get([x as f64 / 20.0, y as f64 / 20.0]) * 10.0;
            heights[x][y] = noise_value as f32;
        }
    }

    return heights;
}

/**
Creates a mesh from a height map, but because we use the separating axis theorem, the height map has to be guaranteed to be a convex shape
So for each little grid of four points we have to create an individual mesh shape
 */
pub fn create_mesh_from_height_map(height_map: Vec<Vec<f32>>, start_pos: Vector2, dx: f32, dz: f32) -> Vec<Box<dyn MeshShape>> {
    let mut all_meshes: Vec<Box<dyn MeshShape>> = Vec::with_capacity((height_map.len() - 1) * (height_map[0].len() - 1));
    for i in 0..height_map.len() - 1 {
        for j in 0..height_map[0].len() - 1 {
            let ground_mesh = GroundMesh::new(
                Vector2::new(start_pos.x + i as f32 * dx, start_pos.y + j as f32 * dz),
                dx,
                dz,
                height_map[i][j],
                height_map[i + 1][j],
                height_map[i + 1][j + 1],
                height_map[i][j + 1],
            );
            all_meshes.push(
                Box::new(ground_mesh) as Box<dyn MeshShape>
            );
        }
    }
    return all_meshes;
}
