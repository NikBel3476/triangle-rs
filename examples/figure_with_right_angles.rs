use plotters::prelude::*;
use triangle_rs::Builder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
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
        .set_switches("pDq")
        .add_polygon(&polygon)
        .build();

    let root = BitMapBackend::new("examples/plots/figure_with_right_angles.png", (800, 600))
        .into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("figure with right angles", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(30f32..37f32, -39f32..-34f32)?;

    chart.configure_mesh().draw()?;

    let line_series = LineSeries::new(
        tri.points
            .chunks(2)
            .map(|point| (point[0] as f32, point[1] as f32)),
        BLUE.stroke_width(2),
    );

    chart
        .draw_series(line_series)?
        .label("points")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));

    tri.triangle_vertex_iter()
        .map(|triangle| vec![triangle[0], triangle[1], triangle[2], triangle[0]])
        .for_each(|triangle| {
            let triangle_series = LineSeries::new(
                triangle
                    .iter()
                    .map(|point| (point.0 as f32, point.1 as f32)),
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

    println!("Figure area = {:#?}", tri.area());

    Ok(())
}
