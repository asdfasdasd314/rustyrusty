# rustyrusty
I want to try out 3D game programming-esque stuff in Rust!

### Description
A project written with Raylib and Rust (btw) for blazingly fast performance. Mostly this was a project for me to learn about how 3D engines and physics work, so I implemented collision detection *from scratch* using the separating axis theorem and its extension into 3 dimensions, and a naive system for terrain generation. There are certainly bugs remaining, I know because I've seen them, and it's nowhere near like a 3D game or a true physics engine or anything, but I'd say the hard stuff, the low level stuff was completed in this little experiment, and this is an experiment I'll gladly put on the shelf.

### Features
- Separating Axis Theorem implementation in 2D and 3D
- Randomly Generated Terrain (not procedural though, but as I said earlier, this is just really easy to do so I'd rather not)
- Collision detection using the separating axis theorem

### What I learned
This is why I care about this project; what did I learn? I learned of course about 3D physics, specifically collision detection, and how that works. I think I'm never going to use that information again, but I think it's really cool that I was able to build my own---reasonably efficient---collision detection system. I understand the design decisions that go into making a game engine, and the minor optimizations you need to do to speed up your code so it doesn't run at 3 FPS. On the meta level, I learned about wrapping types from external libraries, because I went down a long rabbit hole thinking that I needed more precision to make my code work because Raylib ships f32 Vector3s instead of f64s, but that wasn't the issue. And the thing I absolutely hate to think about now: floating point precision. I spent a good two weeks trying to figure the issues in the floating point arithmetic I had to do in my program, and I tried rounding, increasing the precision, caching certain results to reduce floating point calculations, I researched articles, and prompted ChatGPT every way I knew how, and eventually it didn't even come down to that, but I still think it was incredibly useful.

Overall, I would say I am extremely happy with this project, maybe proud as well, but it looks so just not aesthetically pleasing so I can't really say that. I'm so happy that I was able to learn so much and put so many features I didn't understand, but now do, and I had to research certain topics I never would have encountered anywhere else, and I think I really understand the micro complexities of this field a lot better now.
