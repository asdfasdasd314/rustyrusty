use physics::physics::DynamicBody;

use crate::math::float_precision::*;
use crate::physics::*;
use crate::player::*;
use crate::simulation::game::*;
use crate::simulation::render::*;

mod datastructures;
mod math;
mod physics;
mod shapes;
mod simulation;
mod world;

fn main() {
    let player_rect = RectangularPrism::new(Vector3f64::new(0.0, 0.0, 0.0), 1.0, 1.0, 1.0);
    let rb = DynamicBody::new(player_rect.root, Box::new(player_rect));
    let player = Player::new(
        rb.absolute_position,
        10.0,
        0.1,
        CameraType::ThirdPerson(5.0),
        rb,
    );
    let mut game = Game::new(player, vec![], vec![], (640, 480), "Game".to_string());

    game.game_loop();
}
