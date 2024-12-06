use crate::physics::*;
use crate::game::*;
use crate::player::*;
use raylib::prelude::*;

mod game;
mod math_util;
mod physics;
mod player;
mod world_gen;


fn main() {
    let player_rect = RectangularPrism::new(Vector3::new(0.0, 0.0, 0.0), 0.0, 0.0, 0.0);
    //let rect1 = RectangularPrism::new(Vector3::new(2.0, 3.0, 0.0), 10.0, 20.0, 2.0);
    let rect2 = RectangularPrism::new(Vector3::new(5.0, 3.0, 1.0), 12.0, 5.0, 3.0);
    //let solid_body1 = SolidBody::new(Box::new(rect1));
    let solid_body2 = SolidBody::new(Box::new(rect2));
    let rb = SolidBody::new(Box::new(player_rect));
    let player = Player::new(0.01, 0.1, rb);
    let mut game = Game::new(player, vec![solid_body2], (640, 480), "Game".to_string());

    game.game_loop();
}
