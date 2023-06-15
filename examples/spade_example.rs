use plotters::prelude::*;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let polygon = vec![
        Point2::new(54.903146155527025, -3.989636109517439),
        Point2::new(58.19705982671777, -3.989636109517439),
        Point2::new(58.19705982671777, -8.654993391018085),
        Point2::new(60.658311305934646, -8.667238423252995),
        Point2::new(60.658311305934646, -9.377450292877764),
        Point2::new(54.89702363940957, -9.352960228407944),
        Point2::new(54.903146155527025, -3.989636109517439),
    ];

    let mut cdt = ConstrainedDelaunayTriangulation::<Point2<f64>>::new();
    polygon[0..polygon.len() - 1]
        .iter()
        .enumerate()
        .for_each(|(i, point)| {
            let (vertex0, vertex1) = match i >= polygon.len() - 1 {
                true => (
                    Point2::new(point.x, point.y),
                    Point2::new(polygon[0].x, polygon[0].y),
                ),
                false => (
                    Point2::new(point.x, point.y),
                    Point2::new(polygon[i + 1].x, polygon[i + 1].y),
                ),
            };

            cdt.add_constraint_edge(vertex0, vertex1)
                .unwrap_or_else(|err| {
                    panic!("{err}\n{:?}\n{:?}", vertex0, vertex1);
                });
        });

    let root =
        BitMapBackend::new("examples/plots/spade_example.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("spade_example", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(54f32..61f32, -10f32..-3f32)?;

    chart.configure_mesh().draw()?;

    let line_series = LineSeries::new(
        polygon.iter().map(|point| (point.x as f32, point.y as f32)),
        BLUE.stroke_width(2),
    );

    chart
        .draw_series(line_series)?
        .label("points")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));

    cdt.inner_faces()
        .map(|face| {
            let [a, b, c] = face.vertices();
            vec![a, b, c, a]
        })
        .for_each(|triangle| {
            let triangle_series = LineSeries::new(
                triangle
                    .iter()
                    .map(|point| (point.position().x as f32, point.position().y as f32)),
                &RED,
            );
            chart.draw_series(triangle_series).unwrap();
        });

    // tri.triangle_vertex_iter()
    //     .map(|triangle| vec![triangle[0], triangle[1], triangle[2], triangle[0]])
    //     .for_each(|triangle| {
    //         let triangle_series = LineSeries::new(
    //             triangle
    //                 .iter()
    //                 .map(|point| (point.0 as f32, point.1 as f32)),
    //             &RED,
    //         );
    //         chart.draw_series(triangle_series).unwrap();
    //     });

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.0))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    // println!("Figure area = {:#?}", tri.area());

    Ok(())
}
