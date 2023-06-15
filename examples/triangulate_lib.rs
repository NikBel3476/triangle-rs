use plotters::prelude::*;
use triangulate::PolygonList;
use triangulate::{formats, ListFormat};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let polygon: Vec<Vec<[f32; 2]>> = vec![vec![
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

    let mut triangulated_indices = Vec::<[usize; 2]>::new();
    polygon
        .triangulate(formats::IndexedListFormat::new(&mut triangulated_indices).into_fan_format())
        .expect("Triangulation failed");

    println!("{:?}", triangulated_indices);

    // let tri = Builder::new()
    //     .set_switches("p")
    //     .add_polygon(&polygon)
    //     .build();

    let root =
        BitMapBackend::new("examples/plots/triangulate_lib.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("figure with right angles", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(100)
        .y_label_area_size(100)
        .build_cartesian_2d(30f32..37f32, -39f32..-34f32)?;

    chart.configure_mesh().draw()?;

    let line_series = LineSeries::new(
        polygon[0].iter().map(|point| (point[0], point[1])),
        BLUE.stroke_width(2),
    );

    chart
        .draw_series(line_series)?
        .label("points")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(2)));

    triangulated_indices
        .chunks(3)
        .map(|triangle| {
            let vertex0 = polygon.get_vertex(triangle[0]);
            let vertex1 = polygon.get_vertex(triangle[1]);
            let vertex2 = polygon.get_vertex(triangle[2]);
            vec![vertex0, vertex1, vertex2, vertex0]
        })
        .for_each(|triangle| {
            let triangle_series =
                LineSeries::new(triangle.iter().map(|point| (point[0], point[1])), &RED);
            chart.draw_series(triangle_series).unwrap();
        });

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.0))
        .border_style(BLACK)
        .draw()?;

    root.present()?;

    // println!("Figure area = {:#?}", tri.area());

    Ok(())
}
