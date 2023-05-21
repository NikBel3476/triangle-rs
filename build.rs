fn main() {
    println!("cargo:rerun-if-changed=src/triangle.c");

    let files = ["src/triangle.c"];
    let headers_dirs = ["src"];

    let mut builder = cc::Build::new();
    builder
        .files(files.iter())
        .includes(headers_dirs.iter())
        .define("TRILIBRARY", None)
        .warnings(true)
        .extra_warnings(true);

    if cfg!(windows) {
        builder.define("NO_TIMER", None).define("CPU86", None);
    } else {
        builder.define("LINUX", None);
    }

    builder.compile("triangle.a");
}
