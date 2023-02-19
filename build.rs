fn main() {
	println!("cargo:rerun-if-changed=src/triangle.c");

	let files = ["src/triangle.c"];
	let headers_dirs = ["src"];

	cc::Build::new()
		.files(files.iter())
		.includes(headers_dirs.iter())
		.flag("-O")
		.warnings(true)
		.extra_warnings(true)
		.compile("triangle.a");
}
