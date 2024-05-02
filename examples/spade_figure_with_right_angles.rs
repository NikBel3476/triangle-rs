use plotters::prelude::*;
use spade::{ConstrainedDelaunayTriangulation, Point2, Triangulation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let polygon = vec![
        Point2::new(35.97872543334961, -34.659114837646484),
        Point2::new(35.97872543334961, -37.01911163330078),
        Point2::new(33.9708251953125, -37.01911163330078),
        Point2::new(33.9708251953125, -37.219112396240234),
        Point2::new(34.07872772216797, -37.219112396240234),
        Point2::new(34.0787277221679, -38.4352912902832),
        Point2::new(33.15372467041016, -38.4352912902832),
        Point2::new(33.153724670410156, -37.219112396240234),
        Point2::new(33.25210189819336, -37.219112396240234),
        Point2::new(33.25210189819336, -37.01911163330078),
        Point2::new(32.90689468383789, -37.01911163330078),
        Point2::new(32.90689468383789, -37.219112396240234),
        Point2::new(33.003726959228516, -37.219112396240234),
        Point2::new(33.00372695922856, -38.4352912902832),
        Point2::new(32.0787277221679, -38.4352912902832),
        Point2::new(32.07872772216797, -37.219112396240234),
        Point2::new(32.193763732910156, -37.219112396240234),
        Point2::new(32.19376373291015, -37.01911163330078),
        Point2::new(30.50872802734375, -37.01911163330078),
        Point2::new(30.50872802734375, -34.659114837646484),
        Point2::new(35.97872543334961, -34.659114837646484),
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

    let root = BitMapBackend::new(
        "examples/plots/spade_figure_with_right_angles.png",
        (800, 600),
    )
    .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("spade figure with right angles", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(30f32..37f32, -39f32..-34f32)?;

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
