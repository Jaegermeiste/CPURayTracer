#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define _USE_MATH_DEFINES
#include <chrono>
#include <random>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "stb_image.h"
#include "stb_image_write.h"

#ifdef _WIN32
/*
// Poisson looks worse - banding
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::poisson_distribution<int> distribution(5);
*/
inline double drand48()
{
	double randomNum = 0.0;
	do {
		randomNum = ((double)rand() / RAND_MAX);
		//randomNum = static_cast<double>(distribution(generator)) / 10.0;
	} while (randomNum > 1.0);

	return randomNum;
}
#endif // _WIN32

#include "vec3.h"
#include "ray.h"
#include "hitable.h"
#include "hitable_list.h"
#include "material.h"
#include "sphere.h"
#include "camera.h"

typedef unsigned char byte;

#define SCENE_WIDTH		3840
#define SCENE_HEIGHT	2160
#define IMAGE_STRIDE	3
#define NUM_SAMPLES		100
#define MAX_DEPTH		50

/*
double old_hit_sphere(const vec3& center, double radius, const ray& r)
{
	vec3 oc = r.origin() - center;
	double a = dot(r.direction(), r.direction());
	double b = 2.0 * dot(oc, r.direction());
	double c = dot(oc, oc) - radius * radius;
	double discriminant = b * b - 4 * a* c;
	if (discriminant < 0)
	{
		return -1.0;
	}

	return (-b - sqrt(discriminant)) / (2 * a);
}

vec3 old_color(const ray& r)
{
	double t = old_hit_sphere(vec3(0, 0, -1), 0.5, r);
	if (t > 0)
	{
		// Hit
		vec3 N = unit_vector(r.point_at_parameter(t) - vec3(0, 0, -1));
		return 0.5*vec3(N.x() + 1, N.y() + 1, N.z() + 1);
	}

	// Miss
	vec3 unitDirection = unit_vector(r.direction());
	t = 0.5*(unitDirection.y() + 1.0);
	return (1.0 - t)*vec3(1, 1, 1) + t * vec3(0.5, 0.7, 1.0);
}
*/

void dummy()
{
	byte* data = new byte[SCENE_WIDTH*SCENE_HEIGHT * 3]();

	// Generate dummy image
	for (int y = (SCENE_HEIGHT - 1); y >= 0; y--)
	{
		for (int x = 0; x < SCENE_WIDTH; x++)
		{
			vec3 pixel(static_cast<double>(x) / static_cast<double>(SCENE_WIDTH), static_cast<double>(y) / static_cast<double>(SCENE_HEIGHT), 0.2);

			byte ir = static_cast<byte>(255.99f * pixel[0]);
			byte ig = static_cast<byte>(255.99f * pixel[1]);
			byte ib = static_cast<byte>(255.99f * pixel[2]);

			data[(y * SCENE_WIDTH * IMAGE_STRIDE) + (x * IMAGE_STRIDE) + 0] = ir;
			data[(y * SCENE_WIDTH * IMAGE_STRIDE) + (x * IMAGE_STRIDE) + 1] = ig;
			data[(y * SCENE_WIDTH * IMAGE_STRIDE) + (x * IMAGE_STRIDE) + 2] = ib;
		}
	}
	stbi_flip_vertically_on_write(true);
	stbi_write_bmp("outfile.bmp", SCENE_WIDTH, SCENE_HEIGHT, IMAGE_STRIDE /* RGB */, static_cast<void *>(data));
}

vec3 color(const ray& r, hitable* world, int depth)
{
	hit_record record;
	
	if (world->hit(r, 0.001, FLT_MAX, record))
	{
		// Hit
		ray scattered;
		vec3 attenuation;

		if ((depth < MAX_DEPTH) && record.mat_ptr->scatter(r, record, attenuation, scattered))
		{
			// Return color
			return attenuation * color(scattered, world, depth + 1);
		}
		else
		{
			// Return black
			return vec3(0, 0, 0);
		}

		
		vec3 target = record.p + record.normal + random_in_unit_sphere();
		// old return 0.5*vec3(record.normal.x() + 1, record.normal.y() + 1, record.normal.z() + 1);
		return 0.5 *color(ray(record.p, target - record.p), world, depth + 1);
	}

	// Miss
	vec3 unitDirection = unit_vector(r.direction());
	double t = 0.5*(unitDirection.y() + 1.0);
	return (1.0 - t)*vec3(1, 1, 1) + t * vec3(0.5, 0.7, 1.0);
}

hitable *random_scene() {
	int n = 500;
	hitable **list = new hitable*[n + 1];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
	int i = 1;
	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			double choose_mat = drand48();
			vec3 center(a + 0.9*drand48(), 0.2, b + 0.9*drand48());
			if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
				if (choose_mat < 0.8) {  // diffuse
					list[i++] = new sphere(center, 0.2, new lambertian(vec3(drand48()*drand48(), drand48()*drand48(), drand48()*drand48())));
				}
				else if (choose_mat < 0.95) { // metal
					list[i++] = new sphere(center, 0.2,
						new metal(vec3(0.5*(1 + drand48()), 0.5*(1 + drand48()), 0.5*(1 + drand48())), 0.5*drand48()));
				}
				else {  // glass
					list[i++] = new sphere(center, 0.2, new dielectric(1.5));
				}
			}
		}
	}

	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
	list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

	return new hitable_list(list, i);
}

int main()
{
	srand(static_cast<unsigned int>(time(nullptr)));

	byte* data = new byte[SCENE_WIDTH*SCENE_HEIGHT * 3]();
	double aspect = static_cast<double>(SCENE_WIDTH) / static_cast<double>(SCENE_HEIGHT);

	vec3 lookfrom(13, 2, 3);
	vec3 lookat(0, 0, 0);
	vec3 vup(0, 1, 0);


	hitable* list[5];
	double R = cos(M_PI / 4);
	list[0] = new sphere(vec3(-1, 0, -1), 0.5, new dielectric(1.5));
	list[1] = new sphere(vec3(-1, 0, -1), -0.45, new dielectric(1.5));
	list[2] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0, 1, 0)));
	list[3] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.0));
	list[4] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.5, 0.5, 0.5)));

	hitable *world = new hitable_list(list, 5);

	world = random_scene();

	camera cam(lookfrom, lookat, vup, 20, aspect, 0.1, 10.0); // (lookfrom - lookat).length());
	
	// Generate image
	#pragma omp parallel for
	for (int y = (SCENE_HEIGHT - 1); y >= 0; y--)
	{
		for (int x = 0; x < SCENE_WIDTH; x++)
		{
			vec3 pixelColor(0,0,0);

			for (int sample = 0; sample < NUM_SAMPLES; sample++)
			{
				double u = static_cast<double>(x) / static_cast<double>(SCENE_WIDTH);
				double v = static_cast<double>(y) / static_cast<double>(SCENE_HEIGHT);

				ray r = cam.get_ray(u, v);

				pixelColor += color(r, world, 0);
			}
			pixelColor /= static_cast<double>(NUM_SAMPLES);
			pixelColor = vec3(sqrt(pixelColor[0]), sqrt(pixelColor[1]), sqrt(pixelColor[2]));	// Fake gamma correct
			byte ir = static_cast<byte>(255.99f * pixelColor[0]);
			byte ig = static_cast<byte>(255.99f * pixelColor[1]);
			byte ib = static_cast<byte>(255.99f * pixelColor[2]);

			data[(y * SCENE_WIDTH * IMAGE_STRIDE) + (x * IMAGE_STRIDE) + 0] = ir;
			data[(y * SCENE_WIDTH * IMAGE_STRIDE) + (x * IMAGE_STRIDE) + 1] = ig;
			data[(y * SCENE_WIDTH * IMAGE_STRIDE) + (x * IMAGE_STRIDE) + 2] = ib;
		}
	}
	stbi_flip_vertically_on_write(true);
	stbi_write_bmp("outfile_4k.bmp", SCENE_WIDTH, SCENE_HEIGHT, IMAGE_STRIDE /* RGB */, static_cast<void *>(data));

	return 0;
}