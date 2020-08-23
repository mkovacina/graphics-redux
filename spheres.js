
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

const O = {x:0, y:0, z:0};

class Sphere
{
	constructor(center,radius,color)
	{
		this.center = center;
		this.radius = radius;
		this.color = color;
	}
}

function DotProduct(a,b)
{
	return a.x*b.x + a.y*b.y+a.z*b.z;
}

function VectorSubtract(a,b)
{
	return {x:a.x-b.x, y:a.y-b.y, z:a.z-b.z};
}


const s1 = new Sphere({x:0 ,y:-1,z:3}, 1, [255,0  ,0  ]);
const s2 = new Sphere({x:0 ,y:0 ,z:5}, 1, [0  ,0  ,255]);
const s3 = new Sphere({x:-1,y:0 ,z:4}, 1, [0  ,255,0  ]);

Spheres = [s1,s2,s3];

// just a little optimization
// why calculate this for every invocation
// needs more refactoring to be more testable
//  since it relies on global data
const ViewportWidthToCanvasWidthRatio = Vw/Cw;
const ViewportWidthToCanvasHeightRatio = Vh/Ch;

function CanvasToViewport(Cx,Cy)
{
	const Vx = Cx*ViewportWidthToCanvasWidthRatio;
	const Vy = Cy*ViewportWidthToCanvasHeightRatio;
	const Vz = 1;

	return {x:Vx,y:Vy,z:Vz};
}

const id   = ctx.createImageData(1,1);	// only do this once per page
const data = id.data;						// only do this once per page

function DrawPixel(Cx,Cy,color)
{
	// xxx: later optimize with https://stackoverflow.com/questions/4899799/whats-the-best-way-to-set-a-single-pixel-in-an-html5-canvas
	//ctx.fillStyle = 'rgba(0,0,0,0.25)';
	// convert from canvas to screen coordinates
	
	const sx = Cw/2 + Cx;
	const sy = Ch/2 - Cy;
	
	if (sx < 0 || sx >= Cw || sy < 0 || sy >= Ch) return;
 
	var offset = 4*sx + canvas_pitch*sy;
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

function TraceRay(origin, direction, t_min, t_max)
{
	closest_t = Infinity;
	closest_sphere = null;

	// xxx: add a scene object and add spheres to it later
	for( const sphere of Spheres )
	{
		solutions = IntersectRaySphere(origin, direction, sphere);

		t1 = solutions[0];
		t2 = solutions[1];

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
	if (closest_sphere === null)
	{
		return [255,255,255];
	}

	return closest_sphere.color;
}

function loop() 
{
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

	const t0 = performance.now();

	for ( x = -Cw/2; x < Cw/2; x++)
	{
		for( y = -Ch/2; y < Ch/2; y++ )
		{
			D = CanvasToViewport(x,y);
			color = TraceRay(O,D,1,Infinity);
			DrawPixel(x,y,color)
		}
	}

	const t1 = performance.now();
 
	ctx.putImageData(canvas_buffer, 0, 0);

	console.log(`Call to doSomething took ${t1 - t0} milliseconds.`);
	const div = document.getElementById("performance");
	div.innerHTML = `Call to doSomething took ${t1 - t0} milliseconds.`;
}

function MakePoint(x,y,z) { return { x:x, y:y, z:z } }
function test_DotProduct()
{
	const result1 = DotProduct(MakePoint(0,0,0), MakePoint(0,0,0));
	if (result1 != 0) throw "error";
}

function RunTests()
{
	try
	{
		test_DotProduct();
	}
	catch
	{
		const div = document.getElementById("performance");
		div.innerHTML = "error running tests";
		return false;
	}

	return true;
}


const ready = RunTests();

if (ready) loop();


