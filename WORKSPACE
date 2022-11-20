load("@bazel_tools//tools/build_defs/repo:git.bzl", "git_repository")

local_repository(
    name = "tetgen",
    path = "3rdparty/tetgen",
)

local_repository(
    name = "ann",
    path = "3rdparty/ann",
)

local_repository(
    name = "isosurface",
    path = "3rdparty/isosurface_tao",
)

local_repository(
    name = "trimesh",
    path = "3rdparty/trimesh2",
)

local_repository(
    name = "ply",
    path = "3rdparty/ply",
)

git_repository(
    name = "com_github_gflags_gflags",
    remote = "https://github.com/gflags/gflags.git",
    tag = "v2.2.2",
)

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

_RULES_BOOST_COMMIT = "652b21e35e4eeed5579e696da0facbe8dba52b1f"

http_archive(
    name = "com_github_nelhage_rules_boost",
    sha256 = "c1b8b2adc3b4201683cf94dda7eef3fc0f4f4c0ea5caa3ed3feffe07e1fb5b15",
    strip_prefix = "rules_boost-%s" % _RULES_BOOST_COMMIT,
    urls = [
        "https://github.com/nelhage/rules_boost/archive/%s.tar.gz" % _RULES_BOOST_COMMIT,
    ],
)

load("@com_github_nelhage_rules_boost//:boost/boost.bzl", "boost_deps")

boost_deps()
