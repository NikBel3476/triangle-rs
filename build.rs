fn main() {
    println!("cargo:rerun-if-changed=src/triangle.c");

    let files = ["src/triangle.c"];
    let headers_dirs = ["src"];

    let mut build = cc::Build::new();

    if cfg!(unix) {
        build.flag("-DLINUX");
    } else if cfg!(windows) {
        build.flag("-DCPU86");
    }

    build
        .files(files.iter())
        .flag("-DTRILIBRARY")
        .flag("-DNO_TIMER")
        .flag("-DCDT_ONLY")
        .flag("-DREDUCED")
        .flag("-DANSI_DECLARATORS")
        .includes(headers_dirs.iter())
        .opt_level(3)
        .define("REAL", "double")
        .warnings(true)
        .extra_warnings(true)
        .compile("triangle.a");
}
