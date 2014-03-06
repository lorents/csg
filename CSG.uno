using Uno;
using Uno.Collections;
using Uno.Graphics;
using Uno.Scenes;
using Uno.Content;
using Uno.Content.Models;
using Uno.Scenes.Designer;
using Uno.Math;
using Uno.Vector;

// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.
// 
// Example usage:
// 
//     var cube = CSG.cube();
//     var sphere = CSG.sphere({ radius: 1.3 });
//     var polygons = cube.subtract(sphere).toPolygons();
// 
// ## Implementation Details
// 
// All CSG operations are implemented in terms of two functions, `clipTo()` and
// `invert()`, which remove parts of a BSP tree inside another BSP tree and swap
// solid and empty space, respectively. To find the union of `a` and `b`, we
// want to remove everything in `a` inside `b` and everything in `b` inside `a`,
// then combine polygons from `a` and `b` into one solid:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     a.build(b.allPolygons());
// 
// The only tricky part is handling overlapping coplanar polygons in both trees.
// The code above keeps both copies, but we need to keep them in one tree and
// remove them in the other tree. To remove them from `b` we can clip the
// inverse of `b` against `a`. The code for union now looks like this:
// 
//     a.clipTo(b);
//     b.clipTo(a);
//     b.invert();
//     b.clipTo(a);
//     b.invert();
//     a.build(b.allPolygons());
// 
// Subtraction and intersection naturally follow from set operations. If
// union is `A | B`, subtraction is `A - B = ~(~A | B)` and intersection is
// `A & B = ~(~A | ~B)` where `~` is the complement operator.
// 
// ## License
// 
// Original CSG.js Copyright (c) 2011 Evan Wallace (http://madebyevan.com/), under the MIT license.

// # class CSG

// Holds a binary space partition tree representing a 3D solid. Two solids can
// be combined using the `union()`, `subtract()`, and `intersect()` methods.

public class App : Uno.Application
{
	Polygon[] _polygons;
	
	public App()
	{
		var cube = CSG.cube(float3(0,0,0), 20.0f);
		var cube2 = CSG.cube(float3(21,-21,0), 20.0f);
		_polygons = cube.union(cube2).toPolygons();
	}
	
	public override void Draw()
	{
		foreach (var polygon in _polygons)
		{
			for (int i = 0; i+2 < polygon._vertices.Length; i += 3)
			{
				var a = float4(polygon._vertices[i]._pos,1);
				var b = float4(polygon._vertices[i+1]._pos,1);
				var c = float4(polygon._vertices[i+2]._pos,1);
				draw DefaultShading
				{
					float4[] vertices : new float4[3] { a, b, c };
					VertexPosition : vertex_attrib(vertices).XYZ;
					PixelColor: float4(polygon._vertices[i]._normal,1);
					CullFace : PolygonFace.None;
				};
				
			}
		}
	}
}


public class CSG
{
	public Polygon[] _polygons;

	public CSG(Polygon[] polygons)
	{
		if (polygons == null) 
			throw new ArgumentNullException("polygons");
		_polygons = polygons;
	}
	
	// Construct a CSG solid from a list of `CSG.Polygon` instances.
	static CSG fromPolygons(Polygon[] polygons) 
	{
		return new CSG(polygons);
	}

	public CSG clone() 
	{
		var polygons = new Polygon[_polygons.Length];
		for (int i=0; i<_polygons.Length; i++)
			polygons[i] = _polygons[i].clone();
		
		return new CSG(polygons);
	}

	public Polygon[] toPolygons() 
	{
		return _polygons;
	}

	// Return a new CSG solid representing space in either this solid or in the
	// solid `csg`. Neither this solid nor the solid `csg` are modified.
	// 
	//     A.union(B)
	// 
	//     +-------+            +-------+
	//     |       |            |       |
	//     |   A   |            |       |
	//     |    +--+----+   =   |       +----+
	//     +----+--+    |       +----+       |
	//          |   B   |            |       |
	//          |       |            |       |
	//          +-------+            +-------+
	// 
	public CSG union(CSG csg) 
	{
		var a = new Node(clone().toPolygons());
		var b = new Node(csg.clone().toPolygons());
		a.clipTo(b);
		b.clipTo(a);
		b.invert();
		b.clipTo(a);
		b.invert();
		a.build(b.allPolygons());
		return CSG.fromPolygons(a.allPolygons());
	}

	// Return a new CSG solid representing space in this solid but not in the
	// solid `csg`. Neither this solid nor the solid `csg` are modified.
	// 
	//     A.subtract(B)
	// 
	//     +-------+            +-------+
	//     |       |            |       |
	//     |   A   |            |       |
	//     |    +--+----+   =   |    +--+
	//     +----+--+    |       +----+
	//          |   B   |
	//          |       |
	//          +-------+
	// 
	public CSG subtract(CSG csg) 
	{
		var a = new Node(clone().toPolygons());
		var b = new Node(csg.clone().toPolygons());
		a.invert();
		a.clipTo(b);
		b.clipTo(a);
		b.invert();
		b.clipTo(a);
		b.invert();
		a.build(b.allPolygons());
		a.invert();
		return CSG.fromPolygons(a.allPolygons());
	}

	// Return a new CSG solid representing space both this solid and in the
	// solid `csg`. Neither this solid nor the solid `csg` are modified.
	// 
	//     A.intersect(B)
	// 
	//     +-------+
	//     |       |
	//     |   A   |
	//     |    +--+----+   =   +--+
	//     +----+--+    |       +--+
	//          |   B   |
	//          |       |
	//          +-------+
	// 
	public CSG intersect(CSG csg) 
	{
		var a = new Node(clone().toPolygons());
		var b = new Node(csg.clone().toPolygons());
		a.invert();
		b.clipTo(a);
		b.invert();
		a.clipTo(b);
		b.clipTo(a);
		a.build(b.allPolygons());
		a.invert();
		return CSG.fromPolygons(a.allPolygons());
	}

	// Return a new CSG solid with solid and empty space switched. This solid is
	// not modified.
	public CSG inverse() 
	{
		var csg = clone();
		for (int i =0; i<csg._polygons.Length; i++)
			csg._polygons[i].flip();
		return csg;
	}


	// Construct an axis-aligned solid cuboid. Optional parameters are `center` and
	// `radius`, which default to `[0, 0, 0]` and `[1, 1, 1]`. The radius can be
	// specified using a single number or a list of three numbers, one for each axis.
	// 
	// Example code:
	// 
	//     var cube = CSG.cube({
	//       center: [0, 0, 0],
	//       radius: 1
	//     });

	static float3[] vertices = new [] 
	{ 
		float3(-1,-1, -1),
		float3( 1,-1, -1),
		float3( 1, 1, -1),
		float3(-1, 1, -1),
		float3(-1,-1,  1),
		float3( 1,-1,  1),
		float3( 1, 1,  1),
		float3(-1, 1,  1)
	};

	static ushort[] indices = new ushort[]
	{
		0,1,2,2,3,0,
		1,5,6,6,2,1,
		4,7,6,6,5,4,
		0,3,7,7,4,0,
		5,1,0,0,4,5,
		2,6,7,7,3,2,
	};

	public static CSG cube(float3 center, float radius) 
	{
		var r = float3(radius);
		var polygons = new List<Polygon>();
		for (int i =0; i<indices.Length; i+=3)
		{
			var a = (vertices[indices[i+0]]*radius);
			var b = (vertices[indices[i+1]]*radius);
			var c = (vertices[indices[i+2]]*radius);
			var n = Normalize(Cross(a,b));
			var v = new List<Vertex>();
			v.Add(new Vertex(a,n));	
			v.Add(new Vertex(b,n));	
			v.Add(new Vertex(c,n));	
			polygons.Add(new Polygon(v));
		}
		return CSG.fromPolygons(EnumerableExtensions.ToArray(polygons));
	}

/*
// Construct a solid sphere. Optional parameters are `center`, `radius`,
// `slices`, and `stacks`, which default to `[0, 0, 0]`, `1`, `16`, and `8`.
// The `slices` and `stacks` parameters control the tessellation along the
// longitude and latitude directions.
// 
// Example usage:
// 
//     var sphere = CSG.sphere({
//       center: [0, 0, 0],
//       radius: 1,
//       slices: 16,
//       stacks: 8
//     });
CSG.sphere = function(options) {
  options = options || {};
  var c = new CSG.Vector(options.center || [0, 0, 0]);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var stacks = options.stacks || 8;
  var polygons = [], vertices;
  function vertex(theta, phi) {
    theta *= Math.PI * 2;
    phi *= Math.PI;
    var dir = new CSG.Vector(
      Math.cos(theta) * Math.sin(phi),
      Math.cos(phi),
      Math.sin(theta) * Math.sin(phi)
    );
    vertices.push(new CSG.Vertex(c.plus(dir.times(r)), dir));
  }
  for (var i = 0; i < slices; i++) {
    for (var j = 0; j < stacks; j++) {
      vertices = [];
      vertex(i / slices, j / stacks);
      if (j > 0) vertex((i + 1) / slices, j / stacks);
      if (j < stacks - 1) vertex((i + 1) / slices, (j + 1) / stacks);
      vertex(i / slices, (j + 1) / stacks);
      polygons.push(new CSG.Polygon(vertices));
    }
  }
  return CSG.fromPolygons(polygons);
};

// Construct a solid cylinder. Optional parameters are `start`, `end`,
// `radius`, and `slices`, which default to `[0, -1, 0]`, `[0, 1, 0]`, `1`, and
// `16`. The `slices` parameter controls the tessellation.
// 
// Example usage:
// 
//     var cylinder = CSG.cylinder({
//       start: [0, -1, 0],
//       end: [0, 1, 0],
//       radius: 1,
//       slices: 16
//     });
CSG.cylinder = function(options) {
  options = options || {};
  var s = new CSG.Vector(options.start || [0, -1, 0]);
  var e = new CSG.Vector(options.end || [0, 1, 0]);
  var ray = e.minus(s);
  var r = options.radius || 1;
  var slices = options.slices || 16;
  var axisZ = ray.unit(), isY = (Math.abs(axisZ.y) > 0.5);
  var axisX = new CSG.Vector(isY, !isY, 0).cross(axisZ).unit();
  var axisY = axisX.cross(axisZ).unit();
  var start = new CSG.Vertex(s, axisZ.negated());
  var end = new CSG.Vertex(e, axisZ.unit());
  var polygons = [];
  function point(stack, slice, normalBlend) {
    var angle = slice * Math.PI * 2;
    var out = axisX.times(Math.cos(angle)).plus(axisY.times(Math.sin(angle)));
    var pos = s.plus(ray.times(stack)).plus(out.times(r));
    var normal = out.times(1 - Math.abs(normalBlend)).plus(axisZ.times(normalBlend));
    return new CSG.Vertex(pos, normal);
  }
  for (var i = 0; i < slices; i++) {
    var t0 = i / slices, t1 = (i + 1) / slices;
    polygons.push(new CSG.Polygon([start, point(0, t0, -1), point(0, t1, -1)]));
    polygons.push(new CSG.Polygon([point(0, t1, 0), point(0, t0, 0), point(1, t0, 0), point(1, t1, 0)]));
    polygons.push(new CSG.Polygon([end, point(1, t1, 1), point(1, t0, 1)]));
  }
  return CSG.fromPolygons(polygons);
};
	*/

}
	
// # class Vertex

// Represents a vertex of a polygon. Use your own vertex class instead of this
// one to provide additional features like texture coordinates and vertex
// colors. Custom vertex classes need to provide a `pos` property and `clone()`,
// `flip()`, and `interpolate()` methods that behave analogous to the ones
// defined by `CSG.Vertex`. This class provides `normal` so convenience
// functions like `CSG.sphere()` can return a smooth vertex normal, but `normal`
// is not used anywhere else.

public class Vertex
{
	public float3 _pos;
	public float3 _normal;
	
	public Vertex(float3 pos, float3 normal)
	{
		_pos = pos;
		_normal = normal;
	}

	public Vertex clone() 
	{
		return new Vertex(_pos, _normal);
	}

	// Invert all orientation-specific data (e.g. vertex normal). Called when the
	// orientation of a polygon is flipped.
	public void flip() 
	{
		_normal = -_normal;
	}

	// Create a new vertex between this vertex and `other` by linearly
	// interpolating all properties using a parameter of `t`. Subclasses should
	// override this to interpolate additional properties.
	public Vertex Interpolate(Vertex other, float t) 
	{
		return new Vertex(
			Lerp(_pos, other._pos, t),
			Lerp(_normal, other._normal, t));
	}
}

// # class Plane

// Represents a plane in 3D space.

public class Plane
{
	public float3 _normal;
	public float _w;
	
	public Plane(float3 normal, float w) 
	{
		_normal = normal;
		_w = w;
	}

	// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
	// point is on the plane.
	public static double EPSILON = 1e-5;

	public static Plane fromPoints(float3 a, float3 b, float3 c) 
	{
		var n = Normalize(Cross((b-a), (c-a)));
		return new Plane(n, Dot(n,a));
	}

	public Plane clone()
	{
		return new Plane(_normal, _w);
	}

	public void flip() 
	{
		_normal = -_normal;
		_w = -_w;
	}

	// Split `polygon` by this plane if needed, then put the polygon or polygon
	// fragments in the appropriate lists. Coplanar polygons go into either
	// `coplanarFront` or `coplanarBack` depending on their orientation with
	// respect to this plane. Polygons in front or in back of this plane go into
	// either `front` or `back`.
	public void splitPolygon(Polygon polygon, List<Polygon> coplanarFront, List<Polygon> coplanarBack, List<Polygon> front, List<Polygon> back) 
	{
		
		const int COPLANAR = 0;
		const int FRONT = 1;
		const int BACK = 2;
		const int SPANNING = 3;

		// Classify each point as well as the entire polygon into one of the above
		// four classes.
		var polygonType = 0;
		var types = new List<int>();
		for (var i = 0; i < polygon._vertices.Length; i++) 
		{
			var t = Dot(_normal, polygon._vertices[i]._pos) - _w;
			var type = (t < -Plane.EPSILON) 
				? BACK 
				: (t > Plane.EPSILON) 
					? FRONT 
					: COPLANAR;
			polygonType |= type;
			types.Add(type);
		}
		
		// Put the polygon in the correct list, splitting it when necessary.
		switch (polygonType) 
		{
			case COPLANAR:
				(Dot(_normal,polygon._plane._normal) > 0 
					? coplanarFront 
					: coplanarBack).Add(polygon);
				break;
			case FRONT:
				front.Add(polygon);
				break;
			case BACK:
				back.Add(polygon);
				break;
			case SPANNING:
				var f = new List<Vertex>(), b = new List<Vertex>();
				for (var i = 0; i < polygon._vertices.Length; i++) 
				{
					var j = (i + 1) % polygon._vertices.Length;
					var ti = types[i], tj = types[j];
					var vi = polygon._vertices[i], 
						vj = polygon._vertices[j];
					if (ti != BACK) 
						f.Add(vi);
					if (ti != FRONT) 
						b.Add(ti != BACK ? vi.clone() : vi);
					
					if ((ti | tj) == SPANNING) 
					{
						var t = (_w - Dot(_normal, vi._pos)) / Dot(_normal, vj._pos - vi._pos);
						var v = vi.Interpolate(vj, t);
						f.Add(v);
						b.Add(v.clone());
					}
				}
				if (f.Count >= 3) 
					front.Add(new Polygon(f, polygon._shared));

				if (b.Count >= 3) 
					back.Add(new Polygon(b, polygon._shared));
				break;
		}
		
	}

}

// # class Polygon

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
// 
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).

public class Polygon
{
	public Vertex[] _vertices;
	public object _shared;
	public Plane _plane;

	public Polygon(IEnumerable<Vertex> vertices, object shared = null) 
	{
		_vertices = EnumerableExtensions.ToArray(vertices);
		_shared = shared;
		_plane = Plane.fromPoints(_vertices[0]._pos, _vertices[1]._pos, _vertices[2]._pos);
	}

	public Polygon clone() 
	{
		var vertices = new Vertex[_vertices.Length];
		for (int i=0; i<_vertices.Length; i++)
			vertices[i] = _vertices[i].clone();
		return new Polygon(vertices, _shared);
	}

	public void flip() 
	{
		var vertices = new Vertex[_vertices.Length];
		for (int i=0; i<_vertices.Length; i++)
			vertices[i] = _vertices[i];
		
		for (int i=0; i<vertices.Length; i++)
		{
			vertices[i].flip();
			_vertices[_vertices.Length - i - 1] = vertices[i];
		}
		_plane.flip();
	}
}

// # class Node

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.

public class Node 
{
	public Plane _plane;
	public Node _front;
	public Node _back;
	public List<Polygon> _polygons;
	
	public Node() 
		: this (null)
	{
	
	}
	
	public Node(IEnumerable<Polygon> polygons) 
	{
		_plane = null;
		_front = null;
		_back = null;
		_polygons = new List<Polygon>();
		if (polygons != null) build(polygons);
	}

	public Node clone() 
	{
		var node = new Node();
		if (_plane != null)
			node._plane = _plane.clone();
		if (_front != null)
			node._front = _front.clone();
		if (_back != null)
			node._back = _back.clone();
		node._polygons = new List<Polygon>();
		for(int k=0; k<_polygons.Count; k++)
			node._polygons.Add(_polygons[k].clone());
		return node;
	}

	// Convert solid space to empty space and empty space to solid space.
	public void invert() 
	{
		for (var i = 0; i < _polygons.Count; i++) 
		{
			_polygons[i].flip();
		}
		_plane.flip();
		if (_front != null) 
			_front.invert();
		if (_back != null) 
			_back.invert();
		var temp = _front;
		_front = _back;
		_back = temp;
	}

	// Recursively remove all polygons in `polygons` that are inside this BSP tree.
	public List<Polygon> clipPolygons(List<Polygon> polygons) 
	{
		if (_plane == null) 
			return EnumerableExtensions.ToList(polygons);
		
		var front = new List<Polygon>(), 
			back = new List<Polygon>();
		
		foreach (var polygon in polygons)
		{
			_plane.splitPolygon(polygon, front, back, front, back);
		}
		
		if (_front != null) 
			front = _front.clipPolygons(front);
		
		if (_back != null) 
			back = _back.clipPolygons(back);
		else 
			back = new List<Polygon>();
		
		return EnumerableExtensions.ToList(front.Union(back));
	}

	// Remove all polygons in this BSP tree that are inside the other BSP tree `bsp`.
	public void clipTo(Node bsp) 
	{
		_polygons = bsp.clipPolygons(_polygons);
		if (_front != null) 
			_front.clipTo(bsp);
		if (_back != null) 
			_back.clipTo(bsp);
	}

	// Return a list of all polygons in this BSP tree.
	public Polygon[] allPolygons() 
	{
		var polygons = (IEnumerable<Polygon>)_polygons;
		if (_front != null) 
			polygons = EnumerableExtensions.Union(polygons, _front.allPolygons());
		if (_back != null) 
			polygons = EnumerableExtensions.Union(polygons, _back.allPolygons());
		return EnumerableExtensions.ToArray(polygons);
	}

	// Build a BSP tree out of `polygons`. When called on an existing tree, the
	// new polygons are filtered down to the bottom of the tree and become new
	// nodes there. Each set of polygons is partitioned using the first polygon
	// (no heuristic is used to pick a good split).
	public void build(IEnumerable<Polygon> polygons) 
	{
		var count = EnumerableExtensions.Count(polygons) ;
		if (count == 0) return;

		var first = EnumerableExtensions.First(polygons);
		
		if (_plane == null) 
			_plane = first._plane.clone();
		
		var front = new List<Polygon>(), 
			back = new List<Polygon>();

		foreach (var polygon in polygons) 
		{
			_plane.splitPolygon(polygon, _polygons, _polygons, front, back);
		}

		if (front.Count > 0)
		{
			if (_front == null) 
				_front = new Node();
			_front.build(front);
		}
		
		if (back.Count > 0) 
		{
			if (_back == null) 
				_back = new Node();
			_back.build(back);
		}
	}
}
