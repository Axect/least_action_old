use std::ops::{Add, Sub, Mul, Div};
use std::f64::consts::PI;

/// Lagrangian
///
/// # Description
/// `x, v -> L`
pub type Lagrangian1D = fn(f64, f64, &Vec<f64>) -> f64;

#[allow(non_snake_case)]
pub struct Phys1D {
    pub params: Vec<f64>,
    pub lag: Lagrangian1D,
    pub start: i64,
    pub end: i64,
    pub dt: f64,
}

#[allow(non_snake_case)]
impl Phys1D {
    pub fn new(params: Vec<f64>, lag: Lagrangian1D, start: i64, end: i64, dt: f64) -> Self {
        Phys1D {
            params: params,
            lag: lag.to_owned(),
            start,
            end,
            dt
        }
    }

    pub fn L(&self, x: f64, v: f64) -> f64 {
        (self.lag)(x, v, &self.params)
    }

    pub fn calc_one_node(&self, i: i64) -> f64 {
        // v0 = (q1 - q0) / dt
        // x0 = (q1 + q0) / 2
        let v0 = (i - self.start) as f64 / self.dt;
        let x0 = (self.start + i) as f64 / 2f64;

        // v1 = (q2 - q1) / dt
        // x1 = (q2 + q1) / 2
        let v1 = (self.end - i) as f64 / self.dt;
        let x1 = (i + self.end) as f64 / 2f64;

        self.L(x0, v0) + self.L(x1, v1)
    }

    pub fn brute_force_one_node(&self) -> i64 {
        let mut min_val = std::f64::MAX;
        let mut node = 0i64;

        for i in self.start+1 .. self.end+1 {
            let val = self.calc_one_node(i);
            if val < min_val {
                node = i;
                min_val = val;
            }
        }

        node
    }

    pub fn brute_force(&self, m: usize) -> Vec<i64> {
        let mut min_val = std::f64::MAX;
        let mut min_path = vec![0i64; m+2];

        min_path[0] = self.start;
        min_path[m+1] = self.end;

        for q_vec in comb(self.end - self.start - 1, m).into_iter() {
            let mut action = 0f64;

            let q = q_vec.clone().into_iter()
                .map(|x| x + self.start)
                .collect::<Vec<i64>>();

            let v0 = (q[0] - self.start) as f64 / self.dt;
            let x0 = (self.start + q[0]) as f64 / 2f64;
            action += self.L(x0, v0);

            for i in (1 .. m).rev() {
                let vi = (q[i-1] - q[i]) as f64 / self.dt;
                let xi = (q[i] + q[i-1]) as f64 / 2f64;
                action += self.L(xi, vi);
            }
            let qm = q[m-1];
            let vm = (self.end - qm) as f64 / self.dt;
            let xm = (self.end + qm) as f64 / 2f64;
            action += self.L(xm, vm);

            if action < min_val {
                for j in 1 .. (m+1) {
                    min_path[j] = q[j-1];
                }
                min_val = action;
            }
        }

        min_path
    }
//    pub fn brute_force_search(&self) -> (usize, usize, usize) {
//        let top = self.N + 1;
//        let mut min = std::f64::MAX;
//        let mut min_path = (0usize, 0usize, 0usize);
//
//        for i in 1 .. top {
//            for j in (i+1) .. top {
//                for k in (j+1) .. top {
//                    let vk = (top - k) as f64;
//                    let yk = (top + k) as f64 / 2f64;
//
//                    let vj = (k - j) as f64;
//                    let yj = (k + j) as f64 / 2f64;
//
//                    let vi = (j - i) as f64;
//                    let yi = (j + i) as f64 / 2f64;
//
//                    let v0 = i as f64;
//                    let y0 = i as f64 / 2f64;
//
//                    let action = self.L(yk, vk) + self.L(yj, vj) + self.L(yi, vi) + self.L(y0, v0);
//                    if action < min {
//                        min = action;
//                        min_path = (i, j, k);
//                    }
//                }
//            }
//        }
//
//        min_path
//    }
}

#[allow(non_snake_case)]
pub fn divide_and_conquer(p: Phys1D, level: usize) -> Vec<i64> {
    let pivot = p.brute_force_one_node();

    if level == 1 {
        return vec![pivot];
    }

    let p1 = Phys1D {
        params: p.params.clone(),
        lag: p.lag.clone(),
        start: p.start,
        end: pivot,
        dt: p.dt / 2f64,
    };

    let p2 = Phys1D {
        params: p.params.clone(),
        lag: p.lag.clone(),
        start: pivot,
        end: p.end,
        dt: p.dt / 2f64,
    };

    let mut v1 = divide_and_conquer(p1, level/2);
    v1.push(pivot);
    v1.extend(divide_and_conquer(p2, level/2));
    v1
}

pub type Lagrangian2D = fn(Vec2D, Vec2D, (usize, usize), &Vec<f64>) -> f64;

#[derive(Debug, Clone, Copy)]
pub struct Vec2D {
    pub x: f64,
    pub y: f64,
}

impl PartialEq<Vec2D> for Vec2D {
    fn eq(&self, other: &Self) -> bool {
        ((self.x as usize) == (other.x as usize)) && ((self.y as usize) == (other.y as usize))
    }
}

impl Vec2D {
    pub fn new(x: f64, y: f64) -> Self {
        Vec2D { x, y }
    }

    pub fn norm(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2)).sqrt()
    }

    pub fn norm_square(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2)
    }
}

#[allow(non_snake_case)]
pub struct Phys2D {
    N: (usize, usize),
    params: Vec<f64>,
    lag: Lagrangian2D,
    start: Vec2D,
    end: Vec2D
}

#[allow(non_snake_case)]
impl Phys2D {
    pub fn new(n: (usize, usize), params: Vec<f64>, lag: Lagrangian2D, start: Vec2D, end: Vec2D) -> Self {
        Phys2D {
            N: n,
            params,
            lag: lag.to_owned(),
            start,
            end
        }
    }

    pub fn L(&self, x: Vec2D, v: Vec2D) -> f64 {
        (self.lag)(x, v, self.N, &self.params)
    }

    pub fn calc_one_node(&self, k: Vec2D) -> f64 {
        let v0 = &self.start - &k;
        let x0 = (&self.start + &k) / 2f64;

        let v1 = &self.end - &k;
        let x1 = (&self.end + &k) / 2f64;

        self.L(x0, v0) + self.L(x1, v1)
    }

    pub fn brute_force_one_node(&self) -> (usize, usize) {
        let mut min_action = std::f64::MAX;
        let mut node = (0usize, 0usize);

        for i in 0 .. self.N.0 {
            for j in 0 .. self.N.1 {
                if (i,j) == (self.start.x as usize, self.start.y as usize) {
                    continue;
                } else if (i,j) == (self.end.x as usize, self.end.y as usize) {
                    continue;
                } else {
                    let action = self.calc_one_node(Vec2D::new(i as f64, j as f64));
                    if action < min_action {
                        node = (i, j);
                        min_action = action;
                    }
                }
            }
        }
        node
    }

//    pub fn brute_force(&self, m: usize) -> Vec<(usize, usize)> {
//        let q0 = self.start;
//        let qN = self.end;
//
//        let mut min_action = std::f64::MAX;
//        let mut min_path = vec![(0usize, 0usize); m+2];
//
//        min_path[0] = (q0.x as usize, q0.y as usize);
//        min_path[m+1] = (qN.x as usize, qN.y as usize);
//
//        for q_vec in perm((self.N.0+2)*(self.N.1+2), m).into_iter() {
//            let mut action = 0f64;
//
//            let mut q = q_vec.into_iter()
//                .map(|x| x - 1)
//                .map(|x| one_to_two(x, self.N.0))
//                .filter(|v| !(v == self.start) && !(v == self.end))
//                .collect::<Vec<Vec2D>>();
//
//
//        }
//
//        min_path
//    }

    pub fn brute_force_search(&self) -> Vec<(usize, usize)> {
        let q0 = self.start;
        let qN = self.end;

        let mut min_action = std::f64::MAX;
        let mut min_path = vec![(0usize, 0usize); 5];

        min_path[0] = (q0.x as usize, q0.y as usize);
        min_path[4] = (qN.x as usize, qN.y as usize);

        for ix in 0 .. self.N.0+2 {
            for iy in 0 .. self.N.1 + 2 {
                for jx in ix+1 .. self.N.0+2 {
                    for jy in 0 .. self.N.1+2 {
                        for kx in jx+1 .. self.N.0+2 {
                            for ky in 0 .. self.N.1+2 {
                                let i = Vec2D::new(ix as f64, iy as f64);
                                let j = Vec2D::new(jx as f64, jy as f64);
                                let k = Vec2D::new(kx as f64, ky as f64);

                                if i == q0 || i == qN {
                                    continue
                                } else if j == q0 || j == qN {
                                    continue
                                } else if k == q0 || k == qN {
                                    continue
                                } else {
                                    let v0 = &i - &self.start;
                                    let x0 = (&self.start + &i) / 2f64;

                                    let v1 = &j - &i;
                                    let x1 = (&j + &i) / 2f64;

                                    let v2 = &k - &j;
                                    let x2 = (&k + &j) / 2f64;

                                    let v3 = &self.end - &k;
                                    let x3 = (&self.end + &k) / 2f64;

                                    let action = self.L(x0, v0) + self.L(x1, v1) + self.L(x2, v2) + self.L(x3, v3);
                                    if action < min_action {
                                        min_action = action;
                                        min_path[1] = (i.x as usize, i.y as usize);
                                        min_path[2] = (j.x as usize, j.y as usize);
                                        min_path[3] = (k.x as usize, k.y as usize);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        min_path
    }
}


fn comb(n: i64, k: usize) -> Vec<Vec<i64>> {
    let mut p = (1i64 .. (k+1) as i64).collect::<Vec<i64>>();
    let mut result: Vec<Vec<i64>> = Vec::new();

    loop {
        result.push(p.clone());
        let mut i = k;

        while i > 0 && p[i-1] == (n as usize+i-k) as i64 {
            i -= 1;
        }
        if i > 0 {
            p[i-1] += 1;
            for j in i .. k {
                p[j] = p[j-1] + 1;
            }
        } else {
            break;
        }
    }
    result
}

// 2N+4 ..
// N+2  N+3 .. 2N+3
// 0    1   ..  N+1
#[allow(non_snake_case)]
fn one_to_two(n: usize, N: usize) -> Vec2D {
    let x = n / (N+2) + 1;
    let y = n % (N+2);
    Vec2D::new(x as f64, y as f64)
}


impl<'a> Add<&'a Vec2D> for &'a Vec2D {
    type Output = Vec2D;

    fn add(self, rhs: &Vec2D) -> Self::Output {
        Vec2D::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<'a> Sub<&'a Vec2D> for &'a Vec2D {
    type Output = Vec2D;

    fn sub(self, rhs: &Vec2D) -> Self::Output {
        Vec2D::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<'a> Mul<f64> for &'a Vec2D {
    type Output = Vec2D;

    fn mul(self, rhs: f64) -> Self::Output {
        Vec2D::new(self.x * rhs, self.y * rhs)
    }
}

impl<'a> Div<f64> for &'a Vec2D {
    type Output = Vec2D;

    fn div(self, rhs: f64) -> Self::Output {
        Vec2D::new(self.x / rhs, self.y / rhs)
    }
}

impl Add<Vec2D> for Vec2D {
    type Output = Self;

    fn add(self, rhs: Vec2D) -> Self::Output {
        Vec2D::new(self.x+rhs.x, self.y+rhs.y)
    }
}

impl Sub<Vec2D> for Vec2D {
    type Output = Self;

    fn sub(self, rhs: Vec2D) -> Self::Output {
        Vec2D::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl Mul<f64> for Vec2D {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Vec2D::new(self.x * rhs, self.y * rhs)
    }
}

impl Div<f64> for Vec2D {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Vec2D::new(self.x / rhs, self.y / rhs)
    }
}

pub fn dist_square(x: &Vec2D, y: &Vec2D) -> f64 {
    (x - y).norm()
}

pub fn free_body(_y: f64, v: f64, _param: &Vec<f64>) -> f64 {
    0.5 * v.powi(2)
}

pub fn free_fall(y: f64, v: f64, param: &Vec<f64>) -> f64 {
    let g = param[0];
    0.5 * v.powi(2) - g * y
}

pub fn simple_harmonic(x: f64, v: f64, param: &Vec<f64>) -> f64 {
    let k = param[0];
    0.5 * v.powi(2) - 0.5 * k * x.powi(2)
}

pub fn free_body_2d(_r: Vec2D, v: Vec2D, _n: (usize, usize), _param: &Vec<f64>) -> f64 {
    0.5 * v.norm_square()
}

pub fn free_fall_2d(r: Vec2D, v: Vec2D, _n: (usize, usize), param: &Vec<f64>) -> f64 {
    let g = param[0];
    0.5 * v.norm_square() - g * r.y
}

pub fn circular(r: Vec2D, v: Vec2D, _n: (usize, usize), _param: &Vec<f64>) -> f64 {
    0.5 * v.norm_square() + 2000f64*PI/(16f64*r.norm())
}