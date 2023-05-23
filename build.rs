fn main() {
    println!("cargo:rerun-if-changed=src/triangle.c");

    let files = ["src/triangle.c"];
    let headers_dirs = ["src"];

    let mut builder = cc::Build::new();
    builder
        .files(files.iter())
        .includes(headers_dirs.iter())
        .warnings(true)
        .extra_warnings(true);

    if cfg!(windows) {
        builder.define("CPU86", None);
    }
    if cfg!(linux) {
        builder.define("LINUX", None);
    } else {
        builder.define("NO_TIMER", None);
    }

    builder.compile("triangle.a");
}
