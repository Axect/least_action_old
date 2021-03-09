use least_action::*;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();

    let n: i64 = args[2].parse::<i64>().unwrap() * 50;
    let n2: i64 = (n/50) * (2i64.pow(20)-1) ;
    let p = Phys1D::new(
        vec![],
        free_body,
        0,
        n,
        1f64
    );

    let q = Phys1D::new(
        vec![10f64],
        free_fall,
        0,
        n,
        1f64
    );

    let r = Phys1D::new(
        vec![0.1f64],
        simple_harmonic,
        -n/2,
        n/2,
        1f64
    );

    let pd = Phys1D::new(
        vec![],
        free_body,
        0,
        n2,
        4f64
    );

    let qd = Phys1D::new(
        vec![10f64],
        free_fall,
        0,
        n2,
        4f64
    );

    let rd = Phys1D::new(
        vec![0.1f64],
        simple_harmonic,
        -n2/2,
        n2/2,
        4f64
    );

    let p2 = Phys2D::new(
        (19, 19),
        vec![],
        free_body_2d,
        Vec2D::new(0f64, 0f64),
        Vec2D::new(20f64, 20f64),
    );


    let q2 = Phys2D::new(
        (19, 19),
        vec![10f64],
        free_fall_2d,
        Vec2D::new(0f64, 0f64),
        Vec2D::new(20f64, 0f64),
    );

    let r2 = Phys2D::new(
        (19, 19),
        vec![],
        circular,
        Vec2D::new(0f64,20f64),
        Vec2D::new(20f64, 0f64),
    );

    match args[1].parse::<i64>().unwrap() {
        1 => println!("{:?}", p.brute_force(3)),
        2 => println!("{:?}", q.brute_force(3)),
        3 => println!("{:?}", r.brute_force(3)),
        4 => println!("{:?}", divide_and_conquer(pd, 4)),
        5 => println!("{:?}", divide_and_conquer(qd, 4)),
        6 => println!("{:?}", divide_and_conquer(rd, 4)),
        7 => println!("{:?}", p2.brute_force_search()),
        8 => println!("{:?}", q2.brute_force_search()),
        9 => println!("{:?}", r2.brute_force_search()),
        _ => (),
    }

//    println!("{}", p.brute_force_one_node());
//    println!("{}", q.brute_force_one_node());
//    println!("{}", r.brute_force_one_node());
//    println!("{:?}", p.brute_force(3));
//    println!("{:?}", q.brute_force(3));
//    println!("{:?}", r.brute_force(3));

//    let pd = Phys1D {
//        params: vec![],
//        lag: free_body,
//        start: 0,
//        end: 32,
//        dt: 4f64
//    };
//
//    println!("{:?}", divide_and_conquer(pd, 4));


//
//    println!("{:?}", p2.brute_force_search());
//    println!("{:?}", q2.brute_force_search());
//    println!("{:?}", r2.brute_force_search());
}

