use std::collections::HashMap;
use bresenham::Bresenham;
use image::{DynamicImage, GenericImage, GenericImageView};

// All measurements are performed in micro meters (1000000 um == 1m) and Radians
const THREAD_WIDTH: u32 = 600;
const CANVAS_DIAMETER: u32 = 200000;
const NAIL_DIAMETER: u32 = 2000;
const INPUT_IMAGE_PATH: &str = "images/testcat.png";
const NUM_COLOR_CHANNELS: u8 = 3; // Number of toolheads avaliable to perform stitching with
const THREAD_COLORS_AVALIABLE_COUNT: u8 = 1;
const THREAD_COLORS_AVALIABLE: [(u8, u8, u8); THREAD_COLORS_AVALIABLE_COUNT as usize] = [(0, 0, 0)];
const NUM_NAILS: u16 = 120;
const RESOLUTION: u32 = CANVAS_DIAMETER / THREAD_WIDTH;

fn corner_to_cartesian(x: u32, y: u32) -> (i32, i32) {
    return ((x - (RESOLUTION/2)) as i32, (y - (RESOLUTION/2)) as i32 * -1);
}

fn cartesian_to_corner(x: i32, y: i32) -> (u32, u32) {
    return ((x + (RESOLUTION/2) as i32) as u32, ((y - (RESOLUTION as i32/2)) * -1) as u32);
}

// x and y are the true position of the nail
// l and r coords are the position of the thread on the left and right side of the nail (from the perspective of the center of the canvas)
// Stored values are in corner coordinates, not cartesian
#[derive(Eq, PartialEq, Hash)]
struct Nail {
    name: u16,
    xy: (u32, u32),
    lxy: (u32, u32),
    rxy: (u32, u32),
}

// Theta is measured in Radians
impl Nail {
    fn new(theta: f32, name: u16) -> Self {
        let x: i32 = (((RESOLUTION/2)-1) as f32 * theta.cos()) as i32;
        let y: i32 = (((RESOLUTION/2)-1) as f32 * theta.sin()) as i32;

        let nail_edge_pivot: f32 = (1.0-(((NAIL_DIAMETER as f32/2.0).powi(2)))/(2.0*((CANVAS_DIAMETER as f32/2.0).powi(2)))).acos();

        let lx: i32 = (((RESOLUTION/2)-1) as f32 * (theta + nail_edge_pivot).cos()) as i32;
        let ly: i32 = (((RESOLUTION/2)-1) as f32 * (theta + nail_edge_pivot).sin()) as i32;
        let rx: i32 = (((RESOLUTION/2)-1) as f32 * (theta - nail_edge_pivot).cos()) as i32;
        let ry: i32 = (((RESOLUTION/2)-1) as f32 * (theta - nail_edge_pivot).sin()) as i32;

        Self {
            name,
            xy: cartesian_to_corner(x, y),
            lxy: cartesian_to_corner(lx, ly),
            rxy: cartesian_to_corner(rx, ry),
        }
    }
}

#[derive(Eq, PartialEq, Hash, Copy, Clone, Debug)]
enum Handedness {
    Ll,
    Lr,
    Rl,
    Rr,
}

#[derive(Debug)]
struct Stitch {
    pixels: Vec<(u32, u32)>,
    average_color: (u8, u8, u8),
    average_luminance: f32, // 0 to 255
    handedness: Handedness,
}

impl Stitch {
    fn update_averages(&mut self, pixels: &Vec<Vec<Pixel>>) {
        let mut total_color: (u32, u32, u32) = (0, 0, 0);
        for pixel in &self.pixels {
            total_color.0 += pixels[pixel.0 as usize][pixel.1 as usize].rgb.0 as u32;
            total_color.1 += pixels[pixel.0 as usize][pixel.1 as usize].rgb.1 as u32;
            total_color.2 += pixels[pixel.0 as usize][pixel.1 as usize].rgb.2 as u32;
        }
        self.average_color = ((total_color.0 / self.pixels.len() as u32) as u8, (total_color.1 / self.pixels.len() as u32) as u8, (total_color.2 / self.pixels.len() as u32) as u8);

        self.average_luminance = (0.299*(self.average_color.0 as f32).powi(2) + 0.587*(self.average_color.1 as f32).powi(2) + 0.114*(self.average_color.2 as f32).powi(2)).sqrt();
    }
}

struct Pixel {
    rgb: (u8, u8, u8),
    callback_stitches: Vec<(u16, u16, Handedness)>,
}

impl Pixel {
    fn add_stitch(&mut self, nail1: u16, nail2: u16, handedness: Handedness) {
        self.callback_stitches.push((nail1, nail2, handedness));
    }

    fn fire_callbacks(&self, canvas: &mut Canvas) {
        for callback in &self.callback_stitches {
            canvas.update_stitch(callback.0, callback.1, callback.2);
        }
    }
}

struct Canvas {
    nails: Vec<Nail>,
    stitches: HashMap<(u16, u16, Handedness), Stitch>,
    pixels: Vec<Vec<Pixel>>,
}

impl Canvas {
    fn new(target_image: DynamicImage) -> Self {
        let mut nails: Vec<Nail> = Vec::new();
        for i in 0..NUM_NAILS {
            let theta: f32 = (i as f32) * (2.0 * std::f32::consts::PI / (NUM_NAILS as f32));
            nails.push(Nail::new(theta, i));
        }


        let mut pixels: Vec<Vec<Pixel>> = Vec::new();
        for i in 0..RESOLUTION {
            let mut row: Vec<Pixel> = Vec::new();
            for j in 0..RESOLUTION {
                let in_pixel = target_image.get_pixel(i, j);
                if in_pixel[3] != 255 {
                    row.push(Pixel{rgb: (255, 255, 255), callback_stitches: Vec::new()});
                } else {
                    row.push(Pixel{rgb: (in_pixel[0], in_pixel[1], in_pixel[2]), callback_stitches: Vec::new()});
                }
            }
            pixels.push(row);
        }



        let mut stitches: HashMap<(u16, u16, Handedness), Stitch> = HashMap::new();

        for i in 0..NUM_NAILS {
            for j in i+1..NUM_NAILS {
                let nail1: &Nail = &nails[i as usize];
                let nail2: &Nail = &nails[j as usize];

                let mut four_stitches: (Stitch, Stitch, Stitch, Stitch) = (
                    Stitch { pixels: Vec::new(), average_color: (0, 0, 0), average_luminance: 0.0, handedness: Handedness::Ll },
                    Stitch { pixels: Vec::new(), average_color: (0, 0, 0), average_luminance: 0.0, handedness: Handedness::Lr },
                    Stitch { pixels: Vec::new(), average_color: (0, 0, 0), average_luminance: 0.0, handedness: Handedness::Rl },
                    Stitch { pixels: Vec::new(), average_color: (0, 0, 0), average_luminance: 0.0, handedness: Handedness::Rr },
                );
                 

                for cell in Bresenham::new((nail1.lxy.0 as isize, nail1.lxy.1 as isize), (nail2.lxy.0 as isize, nail2.lxy.1 as isize)) {
                    let (x, y) = cell;
                    four_stitches.0.pixels.push((x as u32, y as u32));
                    pixels[x as usize][y as usize].add_stitch(i, j, Handedness::Ll);
                }
                four_stitches.0.update_averages(&pixels);
                stitches.insert((i, j, Handedness::Ll), four_stitches.0);

                for cell in Bresenham::new((nail1.lxy.0 as isize, nail1.lxy.1 as isize), (nail2.rxy.0 as isize, nail2.rxy.1 as isize)) {
                    let (x, y) = cell;
                    four_stitches.1.pixels.push((x as u32, y as u32));
                    pixels[x as usize][y as usize].add_stitch(i, j, Handedness::Lr);
                }
                four_stitches.1.update_averages(&pixels);
                stitches.insert((i, j, Handedness::Lr), four_stitches.1);

                for cell in Bresenham::new((nail1.rxy.0 as isize, nail1.rxy.1 as isize), (nail2.lxy.0 as isize, nail2.lxy.1 as isize)) {
                    let (x, y) = cell;
                    four_stitches.2.pixels.push((x as u32, y as u32));
                    pixels[x as usize][y as usize].add_stitch(i, j, Handedness::Rl);
                }
                four_stitches.2.update_averages(&pixels);
                stitches.insert((i, j, Handedness::Rl), four_stitches.2);

                for cell in Bresenham::new((nail1.rxy.0 as isize, nail1.rxy.1 as isize), (nail2.rxy.0 as isize, nail2.rxy.1 as isize)) {
                    let (x, y) = cell;
                    four_stitches.3.pixels.push((x as u32, y as u32));
                    pixels[x as usize][y as usize].add_stitch(i, j, Handedness::Rr);
                }
                four_stitches.3.update_averages(&pixels);
                stitches.insert((i, j, Handedness::Rr), four_stitches.3);
            }
        }

        Self {
            nails,
            stitches,
            pixels,
        }
    }

    fn get_stitch(&self, nail1: u16, nail2: u16, hand: Handedness) -> &Stitch {
        let key: (u16, u16, Handedness) = if nail1 < nail2 { (nail1, nail2, hand) } else { (nail2, nail1, hand) };
        return self.stitches.get(&key).unwrap();
    }

    fn update_stitch(&mut self, nail1: u16, nail2: u16, hand: Handedness) {
        let key: (u16, u16, Handedness) = if nail1 < nail2 { (nail1, nail2, hand) } else { (nail2, nail1, hand) };
        self.stitches.get_mut(&key).unwrap().update_averages(&self.pixels);
    }
}

fn main() {
    //Load and process image
    let canvas: Canvas = Canvas::new(preprocess_image());
    //Perform algorithm
    algorithm(canvas);
    //Output result
}

fn preprocess_image() -> image::DynamicImage {
    //Open file
    let img: image::DynamicImage = image::open(INPUT_IMAGE_PATH).expect("Failed to open image");

    //Resize image such that 1 pixel == THREAD_WIDTH with a given CANVAS_DIAMETER
    let resized_img: image::DynamicImage = img.resize_exact(RESOLUTION, RESOLUTION, image::imageops::FilterType::Gaussian);

    println!("Image dimensions: {}x{}", resized_img.width(), resized_img.height());

    return resized_img;
}

fn algorithm(mut canvas: Canvas) {
    let mut debug_image = DynamicImage::new_rgb8(RESOLUTION, RESOLUTION);


    // Make debug image white
    for i in 0..canvas.pixels.len() {
        for j in 0..canvas.pixels[i].len() {
            debug_image.put_pixel(i as u32, j as u32, image::Rgba([255,255,255, 255]));
        }
    }

    // for stitch in canvas.stitches.values() {
    //     for pixel in &stitch.pixels {
    //         debug_image.put_pixel(pixel.0, pixel.1, image::Rgba([255,0,0, 255]));
    //     }
    // }

    // for nail in canvas.nails.iter() {
    //     debug_image.put_pixel(nail.xy.0, nail.xy.1, image::Rgba([0,255,0, 255]));
    //     debug_image.put_pixel(nail.lxy.0, nail.lxy.1, image::Rgba([0,0,255, 255]));
    //     debug_image.put_pixel(nail.rxy.0, nail.rxy.1, image::Rgba([0,0,255, 255]));
    // }


    let mut sorted_stitch_keys: Vec<(u16, u16, Handedness)> = canvas.stitches.keys().cloned().collect();
    loop {
        sorted_stitch_keys.sort_by(|a, b| canvas.stitches[b].average_luminance.partial_cmp(&canvas.stitches[a].average_luminance).unwrap());
        let best_stitch_key = sorted_stitch_keys.pop().unwrap();
        let best_stitch = canvas.stitches.get(&best_stitch_key).unwrap();
        
        if best_stitch.average_luminance > 128.0 {
            break;
        }

        for pixel_xy in best_stitch.pixels.clone() {
            let pixel = &mut canvas.pixels[pixel_xy.0 as usize][pixel_xy.1 as usize];
            pixel.rgb = (255, 255, 255);
            for i in 0..pixel.callback_stitches.len() {
                if pixel.callback_stitches.get(i) == Some(&best_stitch_key) {
                    pixel.callback_stitches.remove(i);
                }
            }
            let callbacks: Vec<(u16, u16, Handedness)> = pixel.callback_stitches.clone();
            for callback in callbacks {
                canvas.update_stitch(callback.0, callback.1, callback.2);
            }
            debug_image.put_pixel(pixel_xy.0, pixel_xy.1, image::Rgba([0,0,0, 255]));
        }

    }

    debug_image.save("images/debug.png").unwrap();
}
