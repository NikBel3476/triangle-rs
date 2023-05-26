use std::ffi::c_void;

type REAL = f64;

/* A few forward declarations.                                               */

/* Labels that signify the result of point location.  The result of a        */
/*   search indicates that the point falls in the interior of a triangle, on */
/*   an edge, on a vertex, or outside the mesh.                              */
#[repr(C)]
enum locateresult {
    INTRIANGLE,
    ONEDGE,
    ONVERTEX,
    OUTSIDE,
}
/* Labels that signify the result of vertex insertion.  The result indicates */
/*   that the vertex was inserted with complete success, was inserted but    */
/*   encroaches upon a subsegment, was not inserted because it lies on a     */
/*   segment, or was not inserted because another vertex occupies the same   */
/*   location.                                                               */
#[repr(C)]
enum insertvertexresult {
    SUCCESSFULVERTEX,
    ENCROACHINGVERTEX,
    VIOLATINGVERTEX,
    DUPLICATEVERTEX,
}
/* Labels that signify the result of direction finding.  The result          */
/*   indicates that a segment connecting the two query points falls within   */
/*   the direction triangle, along the left edge of the direction triangle,  */
/*   or along the right edge of the direction triangle.                      */
#[repr(C)]
enum finddirectionresult {
    WITHIN,
    LEFTCOLLINEAR,
    RIGHTCOLLINEAR,
}

/* An oriented triangle:  includes a pointer to a triangle and orientation.  */
/*   The orientation denotes an edge of the triangle.  Hence, there are      */
/*   three possible orientations.  By convention, each edge always points    */
/*   counterclockwise about the corresponding triangle.                      */
#[repr(C)]
struct otri {
    tri: *mut REAL,
    orient: i32, /* Ranges from 0 to 2. */
}

/* An oriented subsegment:  includes a pointer to a subsegment and an        */
/*   orientation.  The orientation denotes a side of the edge.  Hence, there */
/*   are two possible orientations.  By convention, the edge is always       */
/*   directed so that the "side" denoted is the right side of the edge.      */
#[repr(C)]
struct osub {
    ss: *mut REAL,
    ssorient: i32, /* Ranges from 0 to 1. */
}

/* A queue used to store encroached subsegments.  Each subsegment's vertices */
/*   are stored so that we can check whether a subsegment is still the same. */
#[repr(C)]
struct badsubseg {
    encsubseg: *mut *mut REAL, /* An encroached subsegment. */
    subsegorg: *mut REAL,
    subsegdest: *mut REAL, /* Its two vertices. */
}

/* A queue used to store bad triangles.  The key is the square of the cosine */
/*   of the smallest angle of the triangle.  Each triangle's vertices are    */
/*   stored so that one can check whether a triangle is still the same.      */
#[repr(C)]
struct badtriang {
    poortri: *mut *mut REAL, /* A skinny or too-large triangle. */
    key: REAL,               /* cos^2 of smallest (apical) angle. */
    triangorg: *mut REAL,
    triangdest: *mut REAL,
    triangapex: *mut REAL,      /* Its three vertices. */
    nexttriang: *mut badtriang, /* Pointer to next bad triangle. */
}

/* A stack of triangles flipped during the most recent vertex insertion.     */
/*   The stack is used to undo the vertex insertion if the vertex encroaches */
/*   upon a subsegment.                                                      */
#[repr(C)]
struct flipstacker {
    flippedtri: *mut *mut REAL, /* A recently flipped triangle. */
    prevflip: *mut flipstacker, /* Previous flip in the stack. */
}

/* A node in a heap used to store events for the sweepline Delaunay          */
/*   algorithm.  Nodes do not point directly to their parents or children in */
/*   the heap.  Instead, each node knows its position in the heap, and can   */
/*   look up its parent and children in a separate array.  The `eventptr'    */
/*   points either to a `vertex' or to a triangle (in encoded format, so     */
/*   that an orientation is included).  In the latter case, the origin of    */
/*   the oriented triangle is the apex of a "circle event" of the sweepline  */
/*   algorithm.  To distinguish site events from circle events, all circle   */
/*   events are given an invalid (smaller than `xmin') x-coordinate `xkey'.  */
#[repr(C)]
struct event {
    xkey: REAL,
    ykey: REAL,            /* Coordinates of the event. */
    eventptr: *mut c_void, /* Can be a vertex or the location of a circle event. */
    heapposition: i32,     /* Marks this event's position in the heap. */
}

/* A node in the splay tree.  Each node holds an oriented ghost triangle     */
/*   that represents a boundary edge of the growing triangulation.  When a   */
/*   circle event covers two boundary edges with a triangle, so that they    */
/*   are no longer boundary edges, those edges are not immediately deleted   */
/*   from the tree; rather, they are lazily deleted when they are next       */
/*   encountered.  (Since only a random sample of boundary edges are kept    */
/*   in the tree, lazy deletion is faster.)  `keydest' is used to verify     */
/*   that a triangle is still the same as when it entered the splay tree; if */
/*   it has been rotated (due to a circle event), it no longer represents a  */
/*   boundary edge and should be deleted.                                    */
#[repr(C)]
struct splaynode {
    keyedge: otri,      /* Lprev of an edge on the front. */
    keydest: *mut REAL, /* Used to verify that splay node is still live. */
    lchild: *mut splaynode,
    rchild: *mut splaynode, /* Children in splay tree. */
}

/* A type used to allocate memory.  firstblock is the first block of items.  */
/*   nowblock is the block from which items are currently being allocated.   */
/*   nextitem points to the next slab of free memory for an item.            */
/*   deaditemstack is the head of a linked list (stack) of deallocated items */
/*   that can be recycled.  unallocateditems is the number of items that     */
/*   remain to be allocated from nowblock.                                   */
/*                                                                           */
/* Traversal is the process of walking through the entire list of items, and */
/*   is separate from allocation.  Note that a traversal will visit items on */
/*   the "deaditemstack" stack as well as live items.  pathblock points to   */
/*   the block currently being traversed.  pathitem points to the next item  */
/*   to be traversed.  pathitemsleft is the number of items that remain to   */
/*   be traversed in pathblock.                                              */
/*                                                                           */
/* alignbytes determines how new records should be aligned in memory.        */
/*   itembytes is the length of a record in bytes (after rounding up).       */
/*   itemsperblock is the number of items allocated at once in a single      */
/*   block.  itemsfirstblock is the number of items in the first block,      */
/*   which can vary from the others.  items is the number of currently       */
/*   allocated items.  maxitems is the maximum number of items that have     */
/*   been allocated at once; it is the current number of items plus the      */
/*   number of records kept on deaditemstack.                                */
#[repr(C)]
pub struct memorypool {
    pub firstblock: *mut *mut c_void,
    pub nowblock: *mut *mut c_void,
    pub nextitem: *mut c_void,
    pub deaditemstack: *mut c_void,
    pub pathblock: *mut *mut c_void,
    pub pathitem: *mut c_void,
    pub alignbytes: i32,
    pub itembytes: i32,
    pub itemsperblock: i32,
    pub itemsfirstblock: i32,
    pub items: u64,
    pub maxitems: u64,
    pub unallocateditems: i32,
    pub pathitemsleft: i32,
}

// FIXME: not thread-safe(it is neccessary to remove global variables)
/* Global constants.                                                         */
static splitter: REAL = 0.0; /* Used to split REAL factors for exact multiplication. */
static epsilon: REAL = 0.0; /* Floating-point machine epsilon. */
static resulterrbound: REAL = 0.0;
static ccwerrboundA: REAL = 0.0;
static ccwerrboundB: REAL = 0.0;
static ccwerrboundC: REAL = 0.0;
static iccerrboundA: REAL = 0.0;
static iccerrboundB: REAL = 0.0;
static iccerrboundC: REAL = 0.0;
static o3derrboundA: REAL = 0.0;
static o3derrboundB: REAL = 0.0;
static o3derrboundC: REAL = 0.0;

/* Random number seed is not constant, but I've made it global anyway.       */
static randomseed: u64 = 0; /* Current random number seed. */

/* Mesh data structure.  Triangle operates on only one mesh, but the mesh    */
/*   structure is used (instead of global variables) to allow reentrancy.    */
#[repr(C)]
struct mesh {
    /* Variables used to allocate memory for triangles, subsegments, vertices,   */
    /*   viri (triangles being eaten), encroached segments, bad (skinny or too   */
    /*   large) triangles, and splay tree nodes.                                 */
    triangles: memorypool,
    subsegs: memorypool,
    vertices: memorypool,
    viri: memorypool,
    badsubsegs: memorypool,
    badtriangles: memorypool,
    flipstackers: memorypool,
    splaynodes: memorypool,

    /* Variables that maintain the bad triangle queues.  The queues are          */
    /*   ordered from 4095 (highest priority) to 0 (lowest priority).            */
    queuefront: [*mut badtriang; 4096],
    queuetail: [*mut badtriang; 4096],
    nextnonemptyq: [i32; 4096],
    firstnonemptyq: i32,

    /* Variable that maintains the stack of recently flipped triangles.          */
    lastflip: *mut flipstacker,

    /* Other variables. */
    xmin: REAL,
    xmax: REAL,
    ymin: REAL,
    ymax: REAL,           /* x and y bounds. */
    xminextreme: REAL,    /* Nonexistent x value used as a flag in sweepline. */
    invertices: i32,      /* Number of input vertices. */
    inelements: i32,      /* Number of input triangles. */
    insegments: i32,      /* Number of input segments. */
    holes: i32,           /* Number of input holes. */
    regions: i32,         /* Number of input regions. */
    undeads: i32,         /* Number of input vertices that don't appear in the mesh. */
    edges: i64,           /* Number of output edges. */
    mesh_dim: i32,        /* Dimension (ought to be 2). */
    nextras: i32,         /* Number of attributes per vertex. */
    eextras: i32,         /* Number of attributes per triangle. */
    hullsize: i64,        /* Number of edges in convex hull. */
    steinerleft: i32,     /* Number of Steiner points not yet used. */
    vertexmarkindex: i32, /* Index to find boundary marker of a vertex. */
    vertex2triindex: i32, /* Index to find a triangle adjacent to a vertex. */
    highorderindex: i32,  /* Index to find extra nodes for high-order elements. */
    elemattribindex: i32, /* Index to find attributes of a triangle. */
    areaboundindex: i32,  /* Index to find area bound of a triangle. */
    checksegments: i32,   /* Are there segments in the triangulation yet? */
    checkquality: i32,    /* Has quality triangulation begun yet? */
    readnodefile: i32,    /* Has a .node file been read? */
    samples: i64,         /* Number of random samples for point location. */

    incirclecount: i32,     /* Number of incircle tests performed. */
    counterclockcount: i32, /* Number of counterclockwise tests performed. */
    orient3dcount: i32,     /* Number of 3D orientation tests performed. */
    hyperbolacount: i32,    /* Number of right-of-hyperbola tests performed. */
    circumcentercount: i32, /* Number of circumcenter calculations performed. */
    circletopcount: i32,    /* Number of circle top calculations performed. */

    /* Triangular bounding box vertices.                                         */
    infvertex1: *mut REAL,
    infvertex2: *mut REAL,
    infvertex3: *mut REAL,

    /* Pointer to the `triangle' that occupies all of "outer space."             */
    dummytri: *mut REAL,
    dummytribase: *mut REAL, /* Keep base address so we can free() it later. */

    /* Pointer to the omnipresent subsegment.  Referenced by any triangle or     */
    /*   subsegment that isn't really connected to a subsegment at that          */
    /*   location.                                                               */
    dummysub: *mut REAL,
    dummysubbase: *mut REAL, /* Keep base address so we can free() it later. */

    /* Pointer to a recently visited triangle.  Improves point location if       */
    /*   proximate vertices are inserted sequentially.                           */
    recenttri: otri,
} /* End of `struct mesh'. */

/* Data structure for command line switches and file names.  This structure  */
/*   is used (instead of global variables) to allow reentrancy.              */
#[repr(C)]
struct behavior {
    /* Switches for the triangulator.                                            */
    /*   poly: -p switch.  refine: -r switch.                                    */
    /*   quality: -q switch.                                                     */
    /*     minangle: minimum angle bound, specified after -q switch.             */
    /*     goodangle: cosine squared of minangle.                                */
    /*     offconstant: constant used to place off-center Steiner points.        */
    /*   vararea: -a switch without number.                                      */
    /*   fixedarea: -a switch with number.                                       */
    /*     maxarea: maximum area bound, specified after -a switch.               */
    /*   usertest: -u switch.                                                    */
    /*   regionattrib: -A switch.  convex: -c switch.                            */
    /*   weighted: 1 for -w switch, 2 for -W switch.  jettison: -j switch        */
    /*   firstnumber: inverse of -z switch.  All items are numbered starting     */
    /*     from `firstnumber'.                                                   */
    /*   edgesout: -e switch.  voronoi: -v switch.                               */
    /*   neighbors: -n switch.  geomview: -g switch.                             */
    /*   nobound: -B switch.  nopolywritten: -P switch.                          */
    /*   nonodewritten: -N switch.  noelewritten: -E switch.                     */
    /*   noiterationnum: -I switch.  noholes: -O switch.                         */
    /*   noexact: -X switch.                                                     */
    /*   order: element order, specified after -o switch.                        */
    /*   nobisect: count of how often -Y switch is selected.                     */
    /*   steiner: maximum number of Steiner points, specified after -S switch.   */
    /*   incremental: -i switch.  sweepline: -F switch.                          */
    /*   dwyer: inverse of -l switch.                                            */
    /*   splitseg: -s switch.                                                    */
    /*   conformdel: -D switch.  docheck: -C switch.                             */
    /*   quiet: -Q switch.  verbose: count of how often -V switch is selected.   */
    /*   usesegments: -p, -r, -q, or -c switch; determines whether segments are  */
    /*     used at all.                                                          */
    /*                                                                           */
    /* Read the instructions to find out the meaning of these switches.          */
    poly: i32,
    refine: i32,
    quality: i32,
    vararea: i32,
    fixedarea: i32,
    usertest: i32,
    regionattrib: i32,
    convex: i32,
    weighted: i32,
    jettison: i32,
    firstnumber: i32,
    edgesout: i32,
    voronoi: i32,
    neighbors: i32,
    geomview: i32,
    nobound: i32,
    nopolywritten: i32,
    nonodewritten: i32,
    noelewritten: i32,
    noiterationnum: i32,
    noholes: i32,
    noexact: i32,
    conformdel: i32,
    incremental: i32,
    sweepline: i32,
    dwyer: i32,
    splitseg: i32,
    docheck: i32,
    quiet: i32,
    verbose: i32,
    usesegments: i32,
    order: i32,
    nobisect: i32,
    steiner: i32,
    minangle: REAL,
    goodangle: REAL,
    offconstant: REAL,
    maxarea: REAL,
} /* End of `struct behavior'. */

/* Fast lookup arrays to speed some of the mesh manipulation primitives.     */

static plus1mod3: [i32; 3] = [1, 2, 0];
static minus1mod3: [i32; 3] = [2, 0, 1];

/// # int triunsuitable(triorg, tridest, triapex, area)
/// * @param vertex triorg;                              The triangle's origin vertex.
/// * @param vertex tridest;                        The triangle's destination vertex.
/// * @param vertex triapex;                               The triangle's apex vertex.
/// * @param REAL area;                                      The area of the triangle.
/// * @return
#[no_mangle]
pub extern "C" fn triunsuitable(
    triorg: &[REAL],  /* *mut REAL */
    tridest: &[REAL], /* *mut REAL */
    triapex: &[REAL], /* *mut REAL */
    area: f64,
) -> i32 {
    let dxoa = triorg[0] - triapex[0];
    let dyoa = triorg[1] - triapex[1];
    let dxda = tridest[0] - triapex[0];
    let dyda = tridest[1] - triapex[1];
    let dxod = triorg[0] - tridest[0];
    let dyod = triorg[1] - tridest[1];
    /* Find the squares of the lengths of the triangle's three edges. */
    let oalen = dxoa * dxoa + dyoa * dyoa;
    let dalen = dxda * dxda + dyda * dyda;
    let odlen = dxod * dxod + dyod * dyod;
    /* Find the square of the length of the longest edge. */
    let mut maxlen = match dalen > oalen {
        true => dalen,
        false => oalen,
    };
    maxlen = match odlen > maxlen {
        true => odlen,
        false => maxlen,
    };

    match maxlen > 0.05 * (triorg[0] * triorg[0] + triorg[1] * triorg[1]) + 0.02 {
        true => 1,
        false => 0,
    }
}

/*****************************************************************************/
/*                                                                           */
/*  poolzero()   Set all of a pool's fields to zero.                         */
/*                                                                           */
/*  This procedure should never be called on a pool that has any memory      */
/*  allocated to it, as that memory would leak.                              */
/*                                                                           */
/*****************************************************************************/
#[no_mangle]
pub extern "C" fn poolzero(pool: &mut memorypool) {
    pool.firstblock = std::ptr::null_mut();
    pool.nowblock = std::ptr::null_mut();
    pool.nextitem = std::ptr::null_mut();
    pool.deaditemstack = std::ptr::null_mut();
    pool.pathblock = std::ptr::null_mut();
    pool.pathitem = std::ptr::null_mut();
    pool.alignbytes = 0;
    pool.itembytes = 0;
    pool.itemsperblock = 0;
    pool.itemsfirstblock = 0;
    pool.items = 0;
    pool.maxitems = 0;
    pool.unallocateditems = 0;
    pool.pathitemsleft = 0;
}
