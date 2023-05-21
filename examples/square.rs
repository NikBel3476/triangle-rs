use plotters::prelude::*;
use triangle_rs::Builder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let polygon = [-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0];
    let tri = Builder::new()
        .set_switches("pDq")
        .add_polygon(&polygon)
        .build();

    let root = BitMapBackend::new("examples/plots/square.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("square", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(-2f32..2f32, -2f32..2f32)?;

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

    Ok(())
}
