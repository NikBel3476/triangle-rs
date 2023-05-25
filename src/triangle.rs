use std::ffi::c_void;

/* A few forward declarations.                                               */

/* Labels that signify the result of point location.  The result of a        */
/*   search indicates that the point falls in the interior of a triangle, on */
/*   an edge, on a vertex, or outside the mesh.                              */
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
enum finddirectionresult {
    WITHIN,
    LEFTCOLLINEAR,
    RIGHTCOLLINEAR,
}

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

// #[repr(C)]
// struct mesh {
//     /* Variables used to allocate memory for triangles, subsegments, vertices,   */
//     /*   viri (triangles being eaten), encroached segments, bad (skinny or too   */
//     /*   large) triangles, and splay tree nodes.                                 */
//     triangles: memorypool,
//     subsegs: memorypool,
//     vertices: memorypool,
//     viri: memorypool,
//     badsubsegs: memorypool,
//     badtriangles: memorypool,
//     flipstackers: memorypool,
//     splaynodes: memorypool,
//
//     /* Variables that maintain the bad triangle queues.  The queues are          */
//     /*   ordered from 4095 (highest priority) to 0 (lowest priority).            */
//
//     struct badtriang *queuefront[4096];
//     struct badtriang *queuetail[4096];
//     int nextnonemptyq[4096];
//     int firstnonemptyq;
//
//     /* Variable that maintains the stack of recently flipped triangles.          */
//
//     struct flipstacker *lastflip;
//
//     /* Other variables. */
//
//     REAL xmin, xmax, ymin, ymax;                            /* x and y bounds. */
//     REAL xminextreme;      /* Nonexistent x value used as a flag in sweepline. */
//     int invertices;                               /* Number of input vertices. */
//     int inelements;                              /* Number of input triangles. */
//     int insegments;                               /* Number of input segments. */
//     int holes;                                       /* Number of input holes. */
//     int regions;                                   /* Number of input regions. */
//     int undeads;    /* Number of input vertices that don't appear in the mesh. */
//     long long edges;                                     /* Number of output edges. */
//     int mesh_dim;                                /* Dimension (ought to be 2). */
//     int nextras;                           /* Number of attributes per vertex. */
//     int eextras;                         /* Number of attributes per triangle. */
//     long long hullsize;                          /* Number of edges in convex hull. */
//     int steinerleft;                 /* Number of Steiner points not yet used. */
//     int vertexmarkindex;         /* Index to find boundary marker of a vertex. */
//     int vertex2triindex;     /* Index to find a triangle adjacent to a vertex. */
//     int highorderindex;  /* Index to find extra nodes for high-order elements. */
//     int elemattribindex;            /* Index to find attributes of a triangle. */
//     int areaboundindex;             /* Index to find area bound of a triangle. */
//     int checksegments;         /* Are there segments in the triangulation yet? */
//     int checkquality;                  /* Has quality triangulation begun yet? */
//     int readnodefile;                           /* Has a .node file been read? */
//     long long samples;              /* Number of random samples for point location. */
//
//     long long incirclecount;                 /* Number of incircle tests performed. */
//     long long counterclockcount;     /* Number of counterclockwise tests performed. */
//     long long orient3dcount;           /* Number of 3D orientation tests performed. */
//     long long hyperbolacount;      /* Number of right-of-hyperbola tests performed. */
//     long long circumcentercount;  /* Number of circumcenter calculations performed. */
//     long long circletopcount;       /* Number of circle top calculations performed. */
//
//     /* Triangular bounding box vertices.                                         */
//
//     vertex infvertex1, infvertex2, infvertex3;
//
//     /* Pointer to the `triangle' that occupies all of "outer space."             */
//
//     triangle *dummytri;
//     triangle *dummytribase;    /* Keep base address so we can free() it later. */
//
//     /* Pointer to the omnipresent subsegment.  Referenced by any triangle or     */
//     /*   subsegment that isn't really connected to a subsegment at that          */
//     /*   location.                                                               */
//
//     subseg *dummysub;
//     subseg *dummysubbase;      /* Keep base address so we can free() it later. */
//
//     /* Pointer to a recently visited triangle.  Improves point location if       */
//     /*   proximate vertices are inserted sequentially.                           */
//
//     struct otri recenttri;
//
// };

pub extern "C" fn triunsuitable(
    triorg: &[f64],
    tridest: &[f64],
    triapex: &[f64],
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
