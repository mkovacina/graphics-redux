"use strict";
// assumption: fixed viewing position
// assumption: fixed camera orientation
// assumption: camera is located at (0,0,0)
// assumption: on axes and directions
//	- the camera is looking down the positive Z-axis
//	- the y-axis is up
//	- the x-axis is positive to the right
// assumption: 'v' will be used to denote vectors (e.g., vD)
// assumption: "<x,y>' denotes the inner product

// declarations
// Camera
//		- O = (Ox,Oy,Oz)			
//		- camera position (see assumption about position)
//
// Canvas/Pixels
//		- (Cx,Cy) the point on the canvas to draw
//
// Viewport
//		- (Vw,Vh) width and height of the viewport
//		- (Vx,Vy,Vz) the point in the viewport for a pixel on the canvas
//		- viewport width and height
//		- perpendicular to Z-axis
//		- centered on the Z-axis
//		- sides are parallel to the x-axis and y-axis
//
// d
//	- the distance between the camera and the viewport
//							
// Field of View (FOV)
//	- defined by the size of the viewport and the distance (d) from the camera
//	- 60 degree FOV produces decent images
//		Vw=Vh=d=1
//
// P 
//	- any point

const canvas = document.querySelector('canvas');
const ctx = canvas.getContext('2d');

// ???: what is window.innerWidth?
const width = canvas.width;//= window.innerWidth;
const height = canvas.height;// = window.innerHeight;

function MakePoint(x,y,z) { return [ x, y, z ] }

// from the gabrielgambetta.com resource
// implement a double-buffered approach
// remember this from Java...
const canvas_buffer = ctx.getImageData(0, 0, canvas.width, canvas.height);
const canvas_pitch = canvas_buffer.width * 4;

const Cw = width;
const Ch = height;

const Vw = 1;
const Vh = 1;
const  d = 1;

const O = MakePoint(0,0,0);

class Sphere
{
	constructor(center,radius,color,specular)
	{
		this.center = center;
		this.radius = radius;
		this.color = color;
		this.specular = specular;
	}
}

// https://stijndewitt.com/2014/01/26/enums-in-javascript/
const LightType = 
	{
		Unknown: 0,
		Ambient: 1,
		Point: 2,
		Directional: 3,
	};

class Light
{
	type = LightType.Unknown;
	intensity = 0;
	position = null;
	direction = null;

	static CreateAmbientLight(intensity)
	{
		var light = new Light();
		light.type = LightType.Ambient;
		light.intensity = intensity;
		return light;
	}

	static CreatePointLight(intensity,position)
	{
		var light = new Light();
		light.type = LightType.Point;
		light.intensity = intensity;
		light.position = position;
		return light;
	}

	static CreateDirectionalLight(intensity,direction)
	{
		var light = new Light();
		light.type = LightType.Directional;
		light.intensity = intensity;
		light.direction = direction;
		return light;
	}
}

function PointsAreEqual(p1,p2)
{
	"use strict";
	return p1[0] === p2[0] && p1[1] === p2[1] && p1[2] === p2[2];
}

function DotProduct(a,b)
{
	"use strict";
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

function VectorLength(x)
{
	"use strict";
	return Math.sqrt(DotProduct(x,x));
}

function VectorSubtract(a,b)
{
	"use strict";
	return MakePoint(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

function VectorAdd(a,b)
{
	"use strict";
	return MakePoint(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

function VectorMultiplyScalar(v,c)
{
	"use strict";
	return MakePoint(v[0]*c, v[1]*c, v[2]*c);
}

// xxx: replace the color array with a color class??
//      everything isn't "a point"
const s1 = new Sphere(MakePoint(0 ,-1,3), 1,       [255,0  ,0  ], 500 ); // shiny
const s2 = new Sphere(MakePoint(2 ,0 ,4), 1,       [0  ,0  ,255], 500 ); // shiny
const s3 = new Sphere(MakePoint(-2,0 ,4), 1,       [0  ,255,0  ], 10  ); // somewhat shiny
const s4 = new Sphere(MakePoint(0,-5001 ,0), 5000, [255,255,0  ], 1000); // very shiny

const l1 = Light.CreateAmbientLight(0.2);
const l2 = Light.CreatePointLight(0.6, MakePoint(2,1,0));
const l3 = Light.CreateDirectionalLight(0.2, MakePoint(1,4,4));

const Spheres = [s1,s2,s3,s4];
const Lights = [l1,l2,l3];

// just a little optimization
// why calculate this for every invocation
// needs more refactoring to be more testable
//  since it relies on global data
const ViewportWidthToCanvasWidthRatio = Vw/Cw;
const ViewportWidthToCanvasHeightRatio = Vh/Ch;

function CanvasToViewport(Cx,Cy)
{
	"use strict";

	const Vx = Cx*ViewportWidthToCanvasWidthRatio;
	const Vy = Cy*ViewportWidthToCanvasHeightRatio;
	const Vz = 1;

	return MakePoint(Vx,Vy,Vz);
}

const id   = ctx.createImageData(1,1);	// only do this once per page
const data = id.data;						// only do this once per page

function DrawPixel(Cx,Cy,color)
{
	"use strict";

	// xxx: later optimize with https://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
	//ctx.fillStyle = 'rgba(0,0,0,0.25)';
	// convert from canvas to screen coordinates

	const sx = Cw/2 + Cx;
	const sy = Ch/2 - Cy;

	if (sx < 0 || sx >= Cw || sy < 0 || sy >= Ch) return;

	const offset = 4*sx + canvas_pitch*sy;
	canvas_buffer.data[offset+0] = color[0];
	canvas_buffer.data[offset+1] = color[1];
	canvas_buffer.data[offset+2] = color[2];
	canvas_buffer.data[offset+3] = 255; 
	//console.log(`canvas (n${Cx},${Cy}) to screen (${sx},${sy})`);
	//data[0]   = color[0];
	//data[1]   = color[1];
	//data[2]   = color[2];
	//data[3]   = 255;
	//ctx.putImageData( id, sx, sy );     
	//console.log(`canvas (n${Cx},${Cy}) to screen (${sx},${sy}) with color ${id.data}`);
	//ctx.fillRect(sx,sy,1,1);
}

function IntersectRaySphere(origin, direction, sphere)
{
	"use strict";

	const C = sphere.center;
	const r = sphere.radius;
	// O(origin) - C
	const OC = VectorSubtract(origin,C);

	// <direction, direction>
	const k1 = DotProduct(direction,direction);
	// 2*<OC,direction>	
	const k2 = 2*DotProduct(OC,direction);
	// <OC,OC> - r*r
	const k3 = DotProduct(OC,OC)-r*r;

	const discriminant = k2*k2-4*k1*k3;

	if (discriminant < 0)
	{
		return [Infinity,Infinity];
	}

	// minor optimization
	// avoid recomputing
	const k1Times2 = k1 << 1;
	const discriminantSquareRooted = Math.sqrt(discriminant);
	const t1 = (-k2 + discriminantSquareRooted) / k1Times2;
	const t2 = (-k2 - discriminantSquareRooted) / k1Times2;
	return [t1,t2];
}

function ComputeSpecularReflectionComponent(V,N,L,s)
{
	"use strict";

	// V: direction towards the camera (viewing)
	// N: normal vector from the point P
	// L: direction of the light ray
	// s: specular reflection factor

	// i: added intensity
	let i = 0;

	const a = 2*DotProduct(N,L);
	const b = VectorMultiplyScalar(N,a);

	// R: direction of the reflected light
	const R = VectorSubtract(b,L);
	const r_dot_v = DotProduct(R, V);

	if (r_dot_v > 0)
	{
		// adding even more?
		// would have thought that this just got computed for whatever component of
		// intensity was added regardless of the type of light

		//i = light.intensity*Math.pow(r_dot_v/VectorLength(reflectedDirection)*VectorLength(viewingDirection),specularity);
		const lengthR = VectorLength(R);
		const lengthV = VectorLength(V);
		const denominator = lengthR * lengthV;
		// isn't this fun
		// let whatnowA = r_dot_v/lengthR*lengthV;
		// let whatnowB = r_dot_v/(lengthR*lengthV);
		//
		// did you know, that whatnowA != whatnowB??
		// awesome...
		// whatnowB is the correct value AFAIK
		i = Math.pow(r_dot_v/denominator,s);
		//i = Math.pow(r_dot_v/lengthR*lengthV,s);
	}

	return i;
}

function findClosestSphere(spheres,O,D,t_min,t_max)
{
	"use strict";

	// spheres:	the list of spheres in the scene
	//		O:	the origin point to compute the intersection from 
	//		D:	the direction of the ray for look for intersections
	//	t_min:	the minimum distance away an object can be
	//	t_max:	the maximum distance away an object can be
	
	let closest_t = Infinity;
	let closest_sphere = null;

	for( const sphere of spheres )
	{
		const [t1,t2] = IntersectRaySphere(O, D, sphere);

		if ( t1 > t_min && t1 < t_max && t1 < closest_t )
		{
			closest_t = t1;
			closest_sphere = sphere;
		}
		if ( t2 > t_min && t2 < t_max && t2 < closest_t )
		{
			closest_t = t2;
			closest_sphere = sphere;
		}
	}

	return [closest_t,closest_sphere];
}

function ComputeLighting(point, normal, viewingDirection, specularity, lights, spheres)
{
	"use strict";

	let intensity = 0;

	for( const light of lights )
	{
		// ambient lighting is the same at all points
		// an all encompassing light
		// there can be no shadow cast by an ambient light
		// as a shadow is not created nor cast but is a void from light
		if (light.type == LightType.Ambient)
		{
			intensity += light.intensity;
		}
		else
		{
			let rayOfLight;
			let t_max;

			if (light.type == LightType.Point)
			{
				// L = light.position - point
				rayOfLight = VectorSubtract(light.position, point);
				t_max = 1.0;
			}
			else
			{
				// L = light.direction
				rayOfLight = light.direction
				t_max = Infinity;
			}

			// my problem...the findClosestSphere is finding way too many shadows
			//
			// check for shadows
			// xxx: optimization, don't find closest, find any
			const [shadow_t, shadow_sphere] = findClosestSphere(spheres,point,rayOfLight,.001, t_max);

			if (shadow_sphere !== null)
				continue;

			// diffuse light
			const n_dot_l = DotProduct(normal,rayOfLight);

			if (n_dot_l > 0)
			{
				const lengthN = VectorLength(normal);
				const lengthL = VectorLength(rayOfLight);
				const numerator = light.intensity * n_dot_l;
				const denominator = lengthN*lengthL;
				const i = numerator / denominator;
				intensity += i;
			}

			// specular light
			// xxx: what if just less than 0?
			if (specularity != -1 )
			{
				intensity += light.intensity * ComputeSpecularReflectionComponent(viewingDirection, normal, rayOfLight, specularity);
			}
		}
	}

	return intensity;
}

function TraceRay(origin, direction, t_min, t_max)
{
	"use strict";

	// origin: the camera effectively
	// direction: vector, the point we are drawing the ray to

	const [closest_t,closest_sphere] = findClosestSphere(Spheres,origin,direction,t_min, t_max) 

	if (closest_sphere === null)
	{
		return [255,255,255];
	}

	//return closest_sphere.color;

	// compute the intersection
	const P = VectorAdd( origin, VectorMultiplyScalar(direction, closest_t));
	// compute sphere normal at the intersection
	//	normal vector for a point on a sphere is the vector
	//	from the center of the sphere to the point
	//	then we divide by the length to make this a normal vector
	let N = VectorSubtract(P, closest_sphere.center);
	//N = N / VectorLength(N);
	const  x = 1.0 / VectorLength(N);
	N = VectorMultiplyScalar(N, 1.0 / VectorLength(N));

	const Dnegative = VectorMultiplyScalar(direction,-1);

	const intensity = ComputeLighting( P, N, Dnegative, closest_sphere.specular, Lights, Spheres );

	//var color = closest_sphere.color * intensity;
	// no mike, this isn't a real language that you are working with
	// nor is it matlab
	// color is 1x3 vector
	// and javascript will gladly multiply this...
	// ...and give you NaN
	// ....and you didn't write tests for TraceRay yet because......
	const color = VectorMultiplyScalar(closest_sphere.color,intensity);

	return color;
}


function assert( predicate, description )
{
	"use strict";
	if (predicate) return;
	console.assert(predicate, description);
	throw description;
}


function test_DotProduct()
{
	"use strict";

	const result1 = DotProduct(MakePoint(0,0,0), MakePoint(0,0,0));
	if (result1 != 0) throw "error";
}

function test_LightCreation()
{
	"use strict";

	// as needed in the future, adjust the intensity comparisons to account
	// for numerical error
	var l1 = Light.CreateAmbientLight(0.2);
	assert( l1.type == LightType.Ambient, "wrong ambient light type" );
	assert( l1.intensity == 0.2, "wrong ambient intensity");
	assert( l1.position === null, "ambient position not null");
	assert( l1.direction === null, "ambient direction not null");

	var l2 = Light.CreatePointLight(0.2, MakePoint(2,1,0));
	assert( l2.type == LightType.Point, "wrong point light type" );
	assert( l2.intensity == 0.2, "wrong point intensity");
	assert( l2.position !== null, "point position is null");
	assert( PointsAreEqual(l2.position, MakePoint(2,1,0)), "point position not correct");
	assert( l2.direction === null, "point direction not null");

	var l3 = Light.CreateDirectionalLight(0.2, MakePoint(1,2,3));
	assert( l3.type == LightType.Directional, "wrong directional light type" );
	assert( l3.intensity == 0.2, "wrong directional intensity");
	assert( l3.position === null, "directional position not null");
	assert( l3.direction !== null, "directional direction is null");
	assert( PointsAreEqual(l3.direction,MakePoint(1,2,3)), "directional direction not null");
}

function test_VectorLength()
{
	"use strict";
	var vector = MakePoint(0,0,0);
	assert( VectorLength(vector) === 0, "vector (0,0,0) length not zero");	
	var vector = MakePoint(1,0,0);
	assert( VectorLength(vector) === 1, "vector (1,0,0) length not zero");	
	var vector = MakePoint(0,1,0);
	assert( VectorLength(vector) === 1, "vector (1,0,0) length not zero");	
	var vector = MakePoint(0,0,1);
	assert( VectorLength(vector) === 1, "vector (1,0,0) length not zero");	
}

function test_VectorMultiplyScalar()
{
	"use strict";
	var vector = MakePoint(1,2,3);
	var expected = MakePoint(2,4,6);
	var actual = VectorMultiplyScalar(vector,2);

	assert( PointsAreEqual(expected,actual), "didn't scalar multiply correctly");
}

function test_ComputeLighting()
{
	"use strict";

	let delta = 0.000000001
    // basic ambient
    {
    	let light1 = Light.CreateAmbientLight(0.1);
    	let lights = [light1];
    	let p = MakePoint(0,0,0);
    	let n = MakePoint(0,0,0);
    	let v = MakePoint(0,0,0);
    	let s = -1;
    	let actual = ComputeLighting(p,n,v,s,[light1],[]);
    	let expected = 0.1;
    	assert( expected === actual, "failed basic ambient light computation" );
    }
    // basic point 
    {
    	let light1 = Light.CreatePointLight(0.1, MakePoint(1,1,1));
    	let p = MakePoint(0,0,0);
    	let n = MakePoint(1,1,1);
    	let v = MakePoint(0,0,0);
    	let s = -1;
    	let actual = ComputeLighting(p,n,v,s,[light1],[]);
    	let expected = 0.1000;
    	let difference = Math.abs(expected-actual);
    	assert( difference < delta, "failed basic point light computation" );
    }
    // basic directional  
    {
    	let light1 = Light.CreateDirectionalLight(0.1, MakePoint(1,1,1));
    	let p = MakePoint(0,0,0);
    	let n = MakePoint(1,1,1);
    	let v = MakePoint(0,0,0);
    	let s = -1;
    	let actual = ComputeLighting(p,n,v,s,[light1],[]);
    	let expected = 0.1000;
    	let difference = Math.abs(expected-actual);
    	assert( difference  < delta, "failed basic directional light computation" );
    }
    // basic combined, no specular reflection
    {
    	let light1 = Light.CreatePointLight(0.1, MakePoint(1,1,1));
    	let light2 = Light.CreatePointLight(0.1, MakePoint(1,1,1));
    	let light3 = Light.CreateDirectionalLight(0.1, MakePoint(1,1,1));
    	let p = MakePoint(0,0,0);
    	let n = MakePoint(1,1,1);
    	let v = MakePoint(0,0,0);
    	let s = -1;
    	let actual = ComputeLighting(p,n,v,s,[light1,light2,light3],[]);
    	let expected = 0.3;
    	let difference = Math.abs(expected-actual);
    	assert( difference < delta, "failed basic combined light computation" );
    }
    // basic combined, with specular reflection
    {
    	let light1 = Light.CreatePointLight(0.1, MakePoint(1,1,1));
    	let light2 = Light.CreatePointLight(0.1, MakePoint(1,1,1));
    	let light3 = Light.CreateDirectionalLight(0.1, MakePoint(1,1,1));
    	let p = MakePoint(0,0,0);
    	let n = MakePoint(1,1,1);
    	let v = MakePoint(1,1,1);
    	let s = 0;
    	let actual = ComputeLighting(p,n,v,s,[light1,light2,light3],[]);
    	let expected = 0.6;
    	let difference = Math.abs(expected-actual);
    	assert( difference < delta, `failed basic combined light computation with specular expected: ${expected}, actual ${actual} ` );
    }
    // basic directional with sphere
    {
    	let light3 = Light.CreateDirectionalLight(0.1, MakePoint(1,1,1));
    	let p = MakePoint(0,0,0);
    	let n = MakePoint(1,1,1);
    	let v = MakePoint(1,1,1);
    	let s = 1;
		const s1 = new Sphere(MakePoint(.5,.5,.5), 1,       [255,0  ,0  ], 500 ); // shiny
    	let actual = ComputeLighting(p,n,v,s,[light3],[s1]);
    	let expected = 0.0;
    	let difference = Math.abs(expected-actual);
    	assert( difference < delta, "failed basic combined light computation" );
    }
}

function test_ComputeSpecularReflectionComponent()
{
	"use strict";

	let V = MakePoint(1,1,1);
	let N = MakePoint(1,1,1);
	let L = MakePoint(1,1,1);
	let s = 2;

	let actual = ComputeSpecularReflectionComponent(V,N,L,s);
	let expected = 1;

	assert( expected === actual, "specular reflection component not computed correctly");
}

function test_findClosestSphere()
{
	"use strict";

	// so, trying to keep the math simple and keep everything on the [1,1,1] vector.
	// initially had spheres with a radius of 1...so, yeah, while the code worked,
	// it wasn't the inuitive assertion result I was looking for, so shrinking the
	// radius to .1 makes the test easier to read.
	
	const s1 = new Sphere(MakePoint(.50,.50,.50), .1, [255,0  ,0  ], 1 ); // shiny
	const s2 = new Sphere(MakePoint(.25,.25,.25), .1, [0  ,0  ,255], 2 ); // shiny
	const s3 = new Sphere(MakePoint(.75,.75,.75), .1, [0  ,255,0  ], 3  ); // somewhat shiny

	const spheres = [s1,s2,s3];

	const O = MakePoint(0,0,0);
	const D = [1,1,1];
	const t_min = 0.001;
	const t_max = 1;

	const [t,sphere] = findClosestSphere(spheres,O,D,t_min,t_max);

	assert( sphere.specular == 2, "found the wrong sphere");

}

function RunTests()
{
	"use strict";
	try
	{
		test_DotProduct();
		test_VectorLength();
		test_LightCreation();
		test_VectorMultiplyScalar();
		test_ComputeLighting();
		test_ComputeSpecularReflectionComponent();
		test_findClosestSphere();
	}
	catch(error)
	{
		const div = document.getElementById("performance");
		div.innerHTML = "error running tests: "+error;
		return false;
	}

	return true;
}

function loop() 
{
	"use strict";

	// rays start from the camera (O) and go through a section of the viewport
	// setup as a parameter equation w.r.t. 't' to define points on the ray
	// P = O + t(V-O)
	// (V-O) is the direction of the ray, noted vD
	// P = O +t*vD
	//
	// if t < 0, that is behind the camera
	// if t in [0,1], that is between the camera and the viewport
	// if t > 1, that is in the scene
	//
	// to define the points on a sphere, assume sphere with center C and radius r
	// then, if P is a point on the sphere, the distance between the C and P is r.
	// |P-C|=r
	// sqrt(<P-C,P-c>)=r where <x,y> denotes inner product
	// <P-C,P-C>=r^2
	//
	// substitute P from the first equestion into the second
	//

	// xxx: add in a check to assure that all light intensities sum to 1.0
	//		this can be added into the scene class...if you ever create that

	try
	{
		const t0 = performance.now();

		for ( let x = -Cw/2; x < Cw/2; x++)
		{
			for( let y = -Ch/2; y < Ch/2; y++ )
			{
				let D = CanvasToViewport(x,y);
				let color = TraceRay(O,D,1,Infinity);
				assert( !isNaN(color[0]), "color is bad");
				assert( !isNaN(color[1]), "color is bad");
				assert( !isNaN(color[2]), "color is bad");
				DrawPixel(x,y,color)
			}
		}

		const t1 = performance.now();

		ctx.putImageData(canvas_buffer, 0, 0);

		console.log(`Call to doSomething took ${t1 - t0} milliseconds.`);
		const div = document.getElementById("performance");
		div.innerHTML = `Call to doSomething took ${t1 - t0} milliseconds.`;
	}
	catch(error)
	{
		const div = document.getElementById("performance");
		div.innerHTML = "error rendering: "+error;

		throw error;
	}
}


//const ready = RunTests();
const ready = true;

if (ready) loop();

// Tasks
// [ ] create tesst for MakePoint
// [ ] create tests for PointsAreEqual
// [ ] create tests for ComputeLighting
// [ ] create tests for VectorAdd
// [ ] create tests for VectorSubtract
// [ ] generate permuatations of lights and variablesj
// [ ] add a scene object to encapsulate lights, spheres, etc
// [ ] remove any magic numbers
