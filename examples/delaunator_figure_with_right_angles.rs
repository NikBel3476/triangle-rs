use delaunator::triangulate;
use plotters::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let polygon: Vec<Vec<[f64; 2]>> = vec![vec![
        [35.97872543334961, -34.659114837646484],
        [35.97872543334961, -37.01911163330078],
        [33.9708251953125, -37.01911163330078],
        [33.9708251953125, -37.219112396240234],
        [34.07872772216797, -37.219112396240234],
        [34.0787277221679, -38.4352912902832],
        [33.15372467041016, -38.4352912902832],
        [33.153724670410156, -37.219112396240234],
        [33.25210189819336, -37.219112396240234],
        [33.25210189819336, -37.01911163330078],
        [32.90689468383789, -37.01911163330078],
        [32.90689468383789, -37.219112396240234],
        [33.003726959228516, -37.219112396240234],
        [33.00372695922856, -38.4352912902832],
        [32.0787277221679, -38.4352912902832],
        [32.07872772216797, -37.219112396240234],
        [32.193763732910156, -37.219112396240234],
        [32.19376373291015, -37.01911163330078],
        [30.50872802734375, -37.01911163330078],
        [30.50872802734375, -34.659114837646484],
        // [35.97872543334961, -34.659114837646484],
    ]];

    let points = polygon[0].iter().map(|polygon_point| {
        delaunator::Point { x: polygon_point[0], y: polygon_point[1] }
    }).collect::<Vec<delaunator::Point>>();

    let triangulation_result = triangulate(&points);
    let triangulated_indices = triangulation_result.triangles.chunks(3).collect::<Vec<_>>();
    println!("{:?}", triangulated_indices);

    let root =
        BitMapBackend::new("examples/plots/delaunator_figure_with_right_angles.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("delaunator figure with right angles", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(30f32..37f32, -39f32..-34f32)?;

    chart.configure_mesh().draw()?;

    let line_series = LineSeries::new(
        polygon[0]
            .iter()
            .map(|point| (point[0] as f32, point[1] as f32)),
        BLUE.stroke_width(2),
    );

    chart
        .draw_series(line_series)?
        .label("points")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));

    let mut figure_area = 0.0;
    triangulated_indices
        .iter()
        .map(|vertex_indices| {
            let vertex0 = polygon[0][vertex_indices[0]];
            let vertex1 = polygon[0][vertex_indices[1]];
            let vertex2 = polygon[0][vertex_indices[2]];
            vec![vertex0, vertex1, vertex2, vertex0]
        })
        .for_each(|triangle| {
            figure_area += 0.5
                * ((triangle[0][0] - triangle[2][0]) * (triangle[1][1] - triangle[2][1])
                    - (triangle[1][0] - triangle[2][0]) * (triangle[0][1] - triangle[2][1]))
                    .abs();
            let triangle_series = LineSeries::new(
                triangle
                    .iter()
                    .map(|point| (point[0] as f32, point[1] as f32)),
                &RED,
            );
            chart.draw_series(triangle_series).unwrap();
        });

    

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.0))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    println!("Figure area = {:#?}", figure_area);

    Ok(())
}
