//! A Rust library wrapped around the 2D mesh generator and Delaunay triangulator [Triangle](https://www.cs.cmu.edu/~quake/triangle.html)

use rayon::prelude::*;
use serde::Serialize;
use serde_pickle as pkl;
use std::ffi::{c_char, c_int};
use std::{error::Error, ffi::CString, fmt, fs::File, path::Path};

//include!("bindings.rs");

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct triangulateio {
    pub pointlist: *mut f64,
    pub pointattributelist: *mut f64,
    pub pointmarkerlist: *mut c_int,
    pub numberofpoints: c_int,
    pub numberofpointattributes: c_int,
    pub trianglelist: *mut c_int,
    pub triangleattributelist: *mut f64,
    pub trianglearealist: *mut f64,
    pub neighborlist: *mut c_int,
    pub numberoftriangles: c_int,
    pub numberofcorners: c_int,
    pub numberoftriangleattributes: c_int,
    pub segmentlist: *mut c_int,
    pub segmentmarkerlist: *mut c_int,
    pub numberofsegments: c_int,
    pub holelist: *mut f64,
    pub numberofholes: c_int,
    pub regionlist: *mut f64,
    pub numberofregions: c_int,
    pub edgelist: *mut c_int,
    pub edgemarkerlist: *mut c_int,
    pub normlist: *mut f64,
    pub numberofedges: c_int,
}
impl Default for triangulateio {
    fn default() -> Self {
        triangulateio {
            pointlist: std::ptr::null_mut::<f64>(),
            pointattributelist: std::ptr::null_mut::<f64>(),
            pointmarkerlist: std::ptr::null_mut::<i32>(),
            numberofpoints: 0i32,
            numberofpointattributes: 0i32,
            trianglelist: std::ptr::null_mut::<i32>(),
            triangleattributelist: std::ptr::null_mut::<f64>(),
            trianglearealist: std::ptr::null_mut::<f64>(),
            neighborlist: std::ptr::null_mut::<i32>(),
            numberoftriangles: 0i32,
            numberofcorners: 0i32,
            numberoftriangleattributes: 0i32,
            segmentlist: std::ptr::null_mut::<i32>(),
            segmentmarkerlist: std::ptr::null_mut::<i32>(),
            numberofsegments: 0i32,
            holelist: std::ptr::null_mut::<f64>(),
            numberofholes: 0i32,
            regionlist: std::ptr::null_mut::<f64>(),
            numberofregions: 0i32,
            edgelist: std::ptr::null_mut::<i32>(),
            edgemarkerlist: std::ptr::null_mut::<i32>(),
            normlist: std::ptr::null_mut::<f64>(),
            numberofedges: 0i32,
        }
    }
}

extern "C" {
    pub fn triangulate(
        arg1: *mut c_char,
        arg2: *mut triangulateio,
        arg3: *mut triangulateio,
        arg4: *mut triangulateio,
    );
}
extern "C" {
    pub fn trifree(memptr: *mut c_int);
}

/// Delaunay triangulation
#[derive(Default, Debug, Serialize)]
pub struct Delaunay {
    /// Triangulation vertices as [x0,y0,x1,y1,...]
    pub points: Vec<f64>,
    /// Indices in `points.chunks(2)`, the first 3 indices correspond to the vertices of the 1st triangle, the next 3 to the 2nd triangle, etc...
    pub point_markers: Vec<i32>,
    pub triangles: Vec<usize>,
    /// List of triangles neighbors: indices in triangles.chunks(3) (3 integers per triangle)
    pub neighbors: Option<Vec<i32>>,
    /// Edges endpoints: indices in points.chunks(2) (2 integers per edge)
    pub edges: Option<Vec<usize>>,
}

impl Delaunay {
    /// Creates a new empty Delaunay triangulation
    pub fn new() -> Self {
        Self {
            points: vec![],
            point_markers: vec![],
            triangles: vec![],
            neighbors: None,
            edges: None,
        }
    }
    /// Returns the `Delaunay` builder
    pub fn builder() -> Builder {
        Builder::new()
    }
    /// Returns the number of vertices
    pub fn n_vertices(&self) -> usize {
        self.points.len() / 2
    }
    /// Returns the number of Delaunay triangles
    pub fn n_triangles(&self) -> usize {
        self.triangles.len() / 3
    }
    /// Returns an iterator over the vertices, each item is a vertex (x,y) coordinates
    pub fn vertex_iter(&self) -> std::slice::Chunks<'_, f64> {
        self.points.chunks(2)
    }
    /// Returns an iterator over mutable vertices, each item is a vertex (x,y) coordinates
    pub fn vertex_iter_mut(&mut self) -> std::slice::ChunksMut<'_, f64> {
        self.points.chunks_mut(2)
    }
    /// Returns a parallel iterator over the vertices, each item is a vertex (x,y) coordinates
    pub fn vertex_par_iter(&self) -> rayon::slice::Chunks<'_, f64> {
        self.points.par_chunks(2)
    }
    /// Returns an iterator over the triangles, each item is the indices of the vertices in `vertex_iter`
    pub fn triangle_iter(&self) -> std::slice::Chunks<'_, usize> {
        self.triangles.chunks(3)
    }
    /// Returns an iterator over mutable triangles, each item is the indices of the vertices in `vertex_iter`
    pub fn triangle_iter_mut(&mut self) -> std::slice::ChunksMut<'_, usize> {
        self.triangles.chunks_mut(3)
    }
    /// Returns a parallel iterator over the triangles, each item is the indices of the vertices in `vertex_iter`
    pub fn triangle_par_iter(&self) -> rayon::slice::Chunks<'_, usize> {
        self.triangles.par_chunks(3)
    }
    /// Gets node x coordinates
    pub fn x(&self) -> Vec<f64> {
        self.vertex_iter().map(|xy| xy[0]).collect()
    }
    /// Gets node y coordinates
    pub fn y(&self) -> Vec<f64> {
        self.vertex_iter().map(|xy| xy[1]).collect()
    }
    /// Returns an interator over the triangles, each item is a vector of the 3 (x,y) vertices coordinates
    pub fn triangle_vertex_iter(&self) -> impl Iterator<Item = Vec<(f64, f64)>> + '_ {
        let x = self.x();
        let y = self.y();
        self.triangle_iter()
            .map(move |t| t.iter().map(|&i| (x[i], y[i])).collect::<Vec<(f64, f64)>>())
    }

    /// Returns the triangle areas
    pub fn triangle_areas(&self) -> Vec<f64> {
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter()
            .map(|t| self.triangle_area(&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]))
            .collect()
    }

    fn triangle_area(&self, a: &[f64], b: &[f64], c: &[f64]) -> f64 {
        0.5 * ((a[0] - c[0]) * (b[1] - c[1]) - (b[0] - c[0]) * (a[1] - c[1])).abs()
    }

    /// Returns the area covered by the mesh as the sum of the triangle area
    pub fn area(&self) -> f64 {
        self.triangle_areas().iter().sum()
    }
    /// Returns the area covered by the mesh as the sum of the Delaunay triangles area
    pub fn mesh_area(&self) -> f64 {
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            // s + 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs()
            s + self.triangle_area(a, b, c)
        })
    }
    pub fn average(&self, vertices: &[[f64; 3]], data: &[f64]) -> f64 {
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
            let sa = t.iter().fold(0., |m, i| m + data[*i]) / 3 as f64;
            s + ta * sa
        })
    }
    pub fn average_with<F: Fn(f64) -> f64>(
        &self,
        vertices: &[[f64; 3]],
        data: &[f64],
        f: F,
    ) -> f64 {
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
            let sa = t.iter().fold(0., |m, i| m + f(data[*i])) / 3 as f64;
            s + ta * sa
        })
    }
    /// Returns the dot product of two vector `a` and `b` sampled on the Delaunay mesh
    pub fn dot(&self, v_a: &[f64], v_b: &[f64]) -> f64 {
        assert_eq!(v_a.len(), self.n_vertices());
        assert_eq!(v_b.len(), self.n_vertices());
        let vertices: Vec<Vec<f64>> = self.vertex_iter().map(|x| x.to_vec()).collect();
        self.triangle_iter().fold(0., |s, t| {
            let (a, b, c) = (&vertices[t[0]], &vertices[t[1]], &vertices[t[2]]);
            let ta = 0.5 * ((a[0] - c[0]) * (b[1] - a[1]) - (a[0] - b[0]) * (c[1] - a[1])).abs();
            let sa = t.iter().fold(0., |m, i| m + v_a[*i]) / 3 as f64;
            let sb = t.iter().fold(0., |m, i| m + v_b[*i]) / 3 as f64;
            s + ta * sa * sb
        }) / self.area()
    }
    /// Returns true if a point `[x,y]` is inside the triangle given by its index (`triangle_id`) in `triangles_iter`, otherwise returns false
    pub fn is_point_inside_triangle(&self, point: &[f64], triangle_id: usize) -> bool {
        let triangle = self.triangle_iter().nth(triangle_id).unwrap();
        let points: Vec<&[f64]> = self.vertex_iter().collect();

        let d1 = self.sign(point, points[triangle[0]], points[triangle[1]]);
        let d2 = self.sign(point, points[triangle[1]], points[triangle[2]]);
        let d3 = self.sign(point, points[triangle[2]], points[triangle[0]]);
        let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
        let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

        !(has_neg && has_pos)
    }

    pub fn is_point_inside(&self, point: &[f64]) -> bool {
        self.which_contains_point(point).is_some()
    }

    fn sign(&self, p1: &[f64], p2: &[f64], p3: &[f64]) -> f64 {
        (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    }
    /// Finds the index of the triangle in `triangles_iter` that contains the given point `[x,y]`
    pub fn which_contains_point(&self, point: &[f64]) -> Option<usize> {
        (0..self.n_triangles()).find(|&k| self.is_point_inside_triangle(point, k))
    }
    /// Returns the barycentric coordinates of a point `[x,y]` with respect to the triangle that contains it
    ///
    /// The triangle that contains the point is specified with the indices `triangle_ids` in `vertex_iter` of the triangle vertices
    pub fn point_into_barycentric(&self, point: &[f64], triangle_ids: &[usize]) -> [f64; 3] {
        let points: Vec<&[f64]> = self.vertex_iter().collect();
        let v: Vec<&[f64]> = triangle_ids.iter().map(|&i| points[i]).collect();
        let area =
            (v[1][1] - v[2][1]) * (v[0][0] - v[2][0]) + (v[2][0] - v[1][0]) * (v[0][1] - v[2][1]);
        let w0 = ((v[1][1] - v[2][1]) * (point[0] - v[2][0])
            + (v[2][0] - v[1][0]) * (point[1] - v[2][1]))
            / area;
        let w1 = ((v[2][1] - v[0][1]) * (point[0] - v[2][0])
            + (v[0][0] - v[2][0]) * (point[1] - v[2][1]))
            / area;
        [w0, w1, 1. - w0 - w1]
    }
    /// Linearly interpolates at a given point [x,y], values `val_at_vertices` at the Delaunay mesh vertices
    ///
    /// The linear interpolation is based on the barycentric coordinates of the point
    pub fn barycentric_interpolation(&self, point: &[f64], val_at_vertices: &[f64]) -> f64 {
        match self.which_contains_point(point) {
            Some(tid) => {
                let triangle_ids = self.triangle_iter().nth(tid).unwrap();
                let values: Vec<f64> = triangle_ids.iter().map(|&i| val_at_vertices[i]).collect();
                self.point_into_barycentric(point, triangle_ids)
                    .iter()
                    .zip(values.iter())
                    .fold(0., |a, (w, v)| a + w * v)
            }
            None => f64::NAN,
        }
    }
    /*
    pub fn linear_interpolation(&self, point: &[f64], val_at_points: &[f64]) -> f64 {
        match self.contain_point(point) {
            Some(tid) => {
                println!("Triangle #{}", tid);
                let ipts = self.triangles[tid].clone();
                let mut wgts = vec![0f64; 3];
                for i in 0..3 {
                    //let j = if i + 1 == 3 { 0 } else { i + 1 };
                    //let k = if i + 2 == 3 { 0 } else { i + 2 };
                    let j = (i + 1) % 3;
                    let k = (i + 2) % 3;
                    wgts[k] = (self.points[ipts[j]][0] - self.points[ipts[i]][0])
                        * (point[1] - self.points[ipts[i]][1])
                        - (self.points[ipts[j]][1] - self.points[ipts[i]][1])
                            * (point[0] - self.points[ipts[i]][0]);
                }
                println!("weights: {:?}", wgts);
                let sum: f64 = wgts.iter().sum();
                let values: Vec<f64> = ipts.iter().map(|k| val_at_points[*k]).collect();
                println!("values: {:?}", values);
                wgts.iter()
                    .zip(values.iter())
                    .fold(0f64, |a, (w, v)| a + w * v / sum)
            }
            None => std::f64::NAN,
        }
    }
     */
    pub fn dump<T: AsRef<Path>>(&self, filename: T) -> Result<(), Box<dyn Error>> {
        pkl::to_writer(&mut File::create(filename)?, self, true)?;
        Ok(())
    }
}

impl fmt::Display for Delaunay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let areas = self.triangle_areas();
        let areas_min = areas.iter().cloned().fold(f64::INFINITY, f64::min);
        let areas_max = areas.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let areas_sum = areas.iter().sum::<f64>();
        let areas_mean = areas_sum / areas.len() as f64;
        write!(
            f,
            r#"Delaunay triangulation:
 - vertices: {}
 - triangles: {}
 - triangles area: 
   - min : {:.6}
   - max : {:.6}
   - mean: {:.6}
   - sum : {:.6}"#,
            self.n_vertices(),
            self.n_triangles(),
            areas_min,
            areas_max,
            areas_mean,
            areas_sum
        )
    }
}

pub enum TriangulateIO {
    Points(Vec<f64>),
}

/// Delaunay triangulation builder
#[derive(Debug)]
pub struct Builder {
    //triangulate_io: Vec<TriangulateIO>,
    points: Vec<f64>,
    segments: Option<Vec<i32>>,
    n_segments: i32,
    holes: Option<Vec<f64>>,
    n_holes: i32,
    switches: String,
    boundary_marker: i32,
    point_markers: Option<Vec<i32>>,
    segment_markers: Option<Vec<i32>>,
    tri_io: triangulateio,
}
impl Builder {
    /// Creates a new Delaunay triangulation builder
    pub fn new() -> Self {
        Default::default()
    }

    /// Sets the Delaunay mesh `x` and `y` vertices coordinates
    pub fn add_nodes(&mut self, nodes: &[f64]) -> &mut Self {
        self.points.extend(nodes);
        self
    }

    pub fn set_segments(self, x: Vec<i32>, y: Vec<i32>) -> Self {
        assert_eq!(x.len(), y.len(), "x and y are not the same length.");
        //let mut data = self.triangulate_io;
        let n = x.len() as i32;
        let xy = x
            .into_iter()
            .zip(y.into_iter())
            .flat_map(|(x, y)| vec![x, y])
            .collect();
        //data.push(TriangulateIO::Points(xy));
        Self {
            segments: Some(xy),
            n_segments: n,
            ..self
        }
    }

    pub fn add_holes(&mut self, x: f64, y: f64) -> &mut Self {
        match self.holes {
            Some(ref mut h) => {
                h.extend(vec![x, y]);
            }
            None => {
                self.holes = Some(vec![x, y]);
            }
        }
        self
    }

    /// Sets the Delaunay mesh vertices as [x0,y0,x1,y1,...]
    pub fn set_tri_points(self, points: Vec<f64>) -> Self {
        /*let mut data = self.triangulate_io;
        data.push(TriangulateIO::Points(points));*/
        Self { points, ..self }
    }

    /// Adds a closed polygon given its vertices [x1,y1,x2,y2,...]
    pub fn add_polygon(&mut self, vertices: &[f64]) -> &mut Self {
        //let boundary_marker = self.boundary_marker + 1;
        let a = (self.points.len() / 2) as i32;
        let n_segments = (vertices.len() / 2) as i32;
        /*
        let point_markers = match self.point_markers.clone() {
            Some(mut p_m) => {
                p_m.extend(vec![boundary_marker; n_segments as usize]);
                p_m
            }
            None => vec![boundary_marker; n_segments as usize],
        };
        let segment_markers = match self.segment_markers.clone() {
            Some(mut s_m) => {
                s_m.extend(vec![boundary_marker; n_segments as usize]);
                s_m
            }
            None => vec![boundary_marker; n_segments as usize],
        };
        println!("point markers: {:?}", point_markers);*/
        let segments_vertices = (0..n_segments)
            .flat_map(|k| vec![a + k, a + (k + 1) % n_segments])
            .collect::<Vec<i32>>();
        match self.segments {
            Some(ref mut s) => {
                s.extend(segments_vertices);
            }
            None => {
                self.segments = Some(segments_vertices);
            }
        };
        self.points.extend(vertices);
        self
    }

    /// Sets triangulation [switches](https://www.cs.cmu.edu/~quake/triangle.switch.html)
    pub fn set_switches(self, switches: &str) -> Self {
        Self {
            switches: format!("z{}", switches),
            ..self
        }
    }

    /// Compute the Delaunay mesh and returns a `Delaunay` structure
    pub fn build(&mut self) -> Delaunay {
        self.tri_io.numberofpoints = (self.points.len() / 2) as i32;
        self.tri_io.pointlist = self.points.as_mut_ptr();
        if let Some(ref mut s) = self.segments {
            self.tri_io.numberofsegments = (s.len() / 2) as i32;
            self.tri_io.segmentlist = s.as_mut_ptr();
        }
        if let Some(ref mut h) = self.holes {
            self.tri_io.numberofholes = (h.len() / 2) as i32;
            self.tri_io.holelist = h.as_mut_ptr();
        }
        //use TriangulateIO::*;
        let mut delaunay: triangulateio = unsafe { std::mem::zeroed() };
        //        println!("Delaunay triangulation with  switches: {}", self.switches);
        let switches = CString::new(self.switches.as_str()).unwrap();
        unsafe {
            let mut empty_tri: triangulateio = std::mem::zeroed();
            triangulate(
                switches.into_raw(),
                &mut self.tri_io,
                &mut delaunay,
                &mut empty_tri,
            )
        };
        let points: Vec<f64> = unsafe {
            let n = delaunay.numberofpoints as usize * 2;
            std::slice::from_raw_parts(delaunay.pointlist, n).to_vec()
        };
        let point_markers: Vec<i32> = unsafe {
            let n = delaunay.numberofpoints as usize;
            std::slice::from_raw_parts(delaunay.pointmarkerlist, n).to_vec()
        };
        let triangles: Vec<usize> = unsafe {
            let n = delaunay.numberoftriangles as usize * 3;
            std::slice::from_raw_parts(delaunay.trianglelist, n).to_vec()
        }
        .iter()
        .map(|x| *x as usize)
        .collect();
        let neighbors: Option<Vec<i32>> = if self.switches.contains("n") {
            let n = delaunay.numberoftriangles as usize * 3;
            Some(unsafe { std::slice::from_raw_parts(delaunay.neighborlist, n).to_vec() })
        } else {
            None
        };
        let edges: Option<Vec<usize>> = if self.switches.contains("e") {
            let n = delaunay.numberofedges as usize * 2;
            Some(
                unsafe { std::slice::from_raw_parts(delaunay.edgelist, n).to_vec() }
                    .iter()
                    .map(|x| *x as usize)
                    .collect(),
            )
        } else {
            None
        };
        Delaunay {
            points,
            point_markers,
            triangles,
            neighbors,
            edges,
        }
    }
}

impl Default for Builder {
    fn default() -> Self {
        Self {
            //triangulate_io: vec![],
            switches: "z".to_owned(),
            points: vec![],
            segments: None,
            n_segments: 0i32,
            holes: None,
            n_holes: 0i32,
            boundary_marker: 1i32,
            point_markers: None,
            segment_markers: None,
            tri_io: triangulateio::default(),
        }
    }
}

impl From<Vec<f64>> for Builder {
    fn from(points: Vec<f64>) -> Self {
        Self {
            //triangulate_io: vec![TriangulateIO::Points(points)],
            points,
            switches: "z".to_owned(),
            segments: None,
            n_segments: 032,
            holes: None,
            n_holes: 0i32,
            boundary_marker: 1i32,
            point_markers: None,
            segment_markers: None,
            tri_io: triangulateio::default(),
        }
    }
}
impl From<&[f64]> for Builder {
    fn from(points: &[f64]) -> Self {
        Self {
            //triangulate_io: vec![TriangulateIO::Points(points.to_owned())],
            points: points.to_owned(),
            switches: "z".to_owned(),
            segments: None,
            n_segments: 0i32,
            holes: None,
            n_holes: 0i32,
            boundary_marker: 1i32,
            point_markers: None,
            segment_markers: None,
            tri_io: triangulateio::default(),
        }
    }
}

impl fmt::Display for Builder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let nodes_1st_line = format!("{} 2 0 0", self.points.len() / 2);
        let nodes = self
            .points
            .chunks(2)
            .enumerate()
            .map(|(k, xy)| format!("{}  {} {}", k, xy[0], xy[1]))
            .collect::<Vec<String>>()
            .join("\n");
        let segs = match &self.segments {
            Some(s) => {
                let segs_1st_line = format!("{} 0", s.len() / 2);
                let segs = s
                    .chunks(2)
                    .enumerate()
                    .map(|(k, xy)| format!("{}  {} {}", k, xy[0], xy[1]))
                    .collect::<Vec<String>>()
                    .join("\n");
                let holes_1st_line = format!("{}", self.n_holes);
                let holes = match &self.holes {
                    Some(h) => h
                        .chunks(2)
                        .enumerate()
                        .map(|(k, xy)| format!("{}  {} {}", k, xy[0], xy[1]))
                        .collect::<Vec<String>>()
                        .join("\n"),
                    None => "".to_owned(),
                };
                [segs_1st_line, segs, holes_1st_line, holes].join("\n")
            }
            None => "".to_owned(),
        };
        write!(f, "{}", [nodes_1st_line, nodes, segs].join("\n"))
    }
}

/*
impl TriPlot for Delaunay {
    fn mesh<'a, D: DrawingBackend>(
        &self,
        x: &[f64],
        y: &[f64],
        color: [u8; 3],
        chart: &mut ChartContext<'a, D, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    ) -> &Self {
        let color = RGBColor(color[0], color[1], color[2]);
        self.triangle_iter()
            .map(|t| t.iter().map(|&i| (x[i], y[i])).collect::<Vec<(f64, f64)>>())
            .for_each(|v| {
                chart
                    .draw_series(LineSeries::new(
                        v.iter().cycle().take(4).map(|(x, y)| (*x, *y)),
                        &color,
                    ))
                    .unwrap();
            });
        self
    }
    fn map<'a, D: DrawingBackend>(
        &self,
        _x: &[f64],
        _y: &[f64],
        _z: &[f64],
        _chart: &mut ChartContext<'a, D, Cartesian2d<RangedCoordf64, RangedCoordf64>>,
    ) -> &Self {
        self
    }
    fn heatmap(
        &self,
        _x: &[f64],
        _y: &[f64],
        _z: &[f64],
        _range: Range<f64>,
        _config: Option<Config>,
    ) -> Result<(), Box<dyn Error>> {
        Ok(())
    }
}
 */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn area_triangle() {
        let nodes = vec![0.0, -1.0, 1.0, 0.0, 0.0, 1.0];
        let tri = Builder::new().set_switches("Q").add_nodes(&nodes).build();
        assert_eq!(tri.area(), 1.0);
    }

    #[test]
    fn area_square() {
        let tri = Builder::new()
            .set_switches("Q")
            .add_nodes(&[0., 0.])
            .add_polygon(&[1., 0., 0., 1., -1., 0., 0., -1.])
            .build();
        assert_eq!(tri.area(), 2.)
    }

    #[test]
    fn area_parallelogram() {
        let nodes = vec![-2.0, -1.0, 2.0, -1.0, 3.0, 1.0, -1.0, 1.0];
        let tri = Builder::new().set_switches("Q").add_nodes(&nodes).build();
        assert_eq!(tri.area(), 8.0);
    }

    #[test]
    fn area_complex_figure_with_right_angles() {
        let polygon = vec![
            35.97872543334961, // 1
            -34.659114837646484,
            35.97872543334961, // 2
            -37.01911163330078,
            33.9708251953125, // 3
            -37.01911163330078,
            33.9708251953125, // 4
            -37.219112396240234,
            34.07872772216797, // 5
            -37.219112396240234,
            34.0787277221679, // 6
            -38.4352912902832,
            33.15372467041016, // 7
            -38.4352912902832,
            33.153724670410156, // 8
            -37.219112396240234,
            33.25210189819336, // 9
            -37.219112396240234,
            33.25210189819336, // 10
            -37.01911163330078,
            32.90689468383789, // 11
            -37.01911163330078,
            32.90689468383789, // 12
            -37.219112396240234,
            33.003726959228516, // 13
            -37.219112396240234,
            33.00372695922856, // 14
            -38.4352912902832,
            32.0787277221679, // 15
            -38.4352912902832,
            32.07872772216797, // 16
            -37.219112396240234,
            32.193763732910156, // 17
            -37.219112396240234,
            32.19376373291015, // 18
            -37.01911163330078,
            30.50872802734375, // 19
            -37.01911163330078,
            30.50872802734375, // 20
            -34.659114837646484,
            35.97872543334961, // 21
            -34.659114837646484,
        ];
        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();
        println!("{:#?}", tri);
        assert_eq!(tri.area(), 15.44548203003071);
    }

    #[test]
    fn point_inside_triangle() {
        let polygon = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();
        let points = vec![
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.5, 0.5],
            [0.0, 1.0],
            [0.0, 0.5],
        ];
        for triangle_index in 0..tri.n_triangles() {
            for point in &points {
                assert!(
                    tri.is_point_inside_triangle(point, triangle_index),
                    "\nPoint {:?} should be inside triangle\n",
                    point
                )
            }
        }
    }

    #[test]
    fn point_outside_triangle() {
        let polygon = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();
        let points = vec![
            [-1.0, -1.0],
            [0.5, -1.0],
            [1.5, -0.5],
            [1.0, 1.0],
            [-0.5, 1.5],
            [-0.5, 0.5],
        ];
        for triangle_index in 0..tri.n_triangles() {
            for point in &points {
                assert!(
                    !tri.is_point_inside_triangle(point, triangle_index),
                    "\nPoint {:?} should be outside triangle\n",
                    point
                )
            }
        }
    }

    #[test]
    fn point_inside_square() {
        let polygon = vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0];
        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();
        let points = vec![
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [1.0, 0.5],
            [1.0, 1.0],
            [0.5, 1.0],
            [0.0, 1.0],
            [0.0, 0.5],
            [0.5, 0.5],
        ];
        for point in &points {
            assert!(
                tri.is_point_inside(point),
                "\nPoint {:?} should be inside square\n",
                point
            )
        }
    }

    #[test]
    fn point_outside_square() {
        let polygon = vec![0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0];
        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();
        let points = vec![
            [-0.5, -0.5],
            [0.5, -0.5],
            [1.5, -0.5],
            [1.5, 0.5],
            [1.5, 1.5],
            [0.5, 1.5],
            [-0.5, 1.5],
            [-0.5, 0.5],
        ];
        for point in &points {
            assert!(
                !tri.is_point_inside(point),
                "\nPoint {:?} should be outside square\n",
                point
            )
        }
    }

    #[test]
    fn rectangles_intersection() {
        let polygon = vec![
            1.0804720876503242,
            9.784116583159095,
            9.452596550210751,
            9.830117267019318,
            9.475596892140864,
            1.2969904109481103,
            1.1034724295804352,
            1.2969904109481103,
            1.0804720876503242,
            9.784116583159095,
        ];

        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();

        let points_outside = vec![
            [7.198563041059868, 10.888132995804426],
            [8.877588001957976, 10.934133679664647],
        ];
        let points_inside = vec![
            [8.854587660027866, 9.577113505788097],
            [7.198563041059868, 9.554113163857984],
        ];

        points_outside.iter().for_each(|point| {
            assert!(
                !tri.is_point_inside(point),
                "\nPoint {:?} should be outside rectangle\n",
                point
            )
        });

        points_inside.iter().for_each(|point| {
            assert!(
                tri.is_point_inside(point),
                "\nPoint {:?} should be inside rectangle\n",
                point
            )
        });
    }

    #[test]
    fn intersection_of_a_figure_and_a_rectangle() {
        let polygon = vec![
            35.97872543334961, // 1
            -34.659114837646484,
            35.97872543334961, // 2
            -37.01911163330078,
            33.9708251953125, // 3
            -37.01911163330078,
            33.9708251953125, // 4
            -37.219112396240234,
            34.07872772216797, // 5
            -37.219112396240234,
            34.0787277221679, // 6
            -38.4352912902832,
            33.15372467041016, // 7
            -38.4352912902832,
            33.153724670410156, // 8
            -37.219112396240234,
            33.25210189819336, // 9
            -37.219112396240234,
            33.25210189819336, // 10
            -37.01911163330078,
            32.90689468383789, // 11
            -37.01911163330078,
            32.90689468383789, // 12
            -37.219112396240234,
            33.003726959228516, // 13
            -37.219112396240234,
            33.00372695922856, // 14
            -38.4352912902832,
            32.0787277221679, // 15
            -38.4352912902832,
            32.07872772216797, // 16
            -37.219112396240234,
            32.193763732910156, // 17
            -37.219112396240234,
            32.19376373291015, // 18
            -37.01911163330078,
            30.50872802734375, // 19
            -37.01911163330078,
            30.50872802734375, // 20
            -34.659114837646484,
            35.97872543334961, // 21
            -34.659114837646484,
        ];

        let tri = Builder::new()
            .set_switches("pDqQ")
            .add_polygon(&polygon)
            .build();

        let points_outside = vec![
            [31.87872886657715, -38.24702072143555],
            [31.87872886657715, -37.34701919555664],
        ];
        let points_inside = vec![
            [32.07872772216797, -38.24702072143555],
            [32.07872772216797, -37.34701919555664],
        ];

        points_outside.iter().for_each(|point| {
            assert!(
                !tri.is_point_inside(point),
                "\nPoint {:?} should be outside figure\n",
                point
            )
        });

        points_inside.iter().for_each(|point| {
            assert!(
                tri.is_point_inside(point),
                "\nPoint {:?} should be inside figure\n",
                point
            )
        });
    }
}
