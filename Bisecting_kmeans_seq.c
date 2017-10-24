#include<stdio.h>
#include<stdlib.h>
#include<math.h>

////Amitsingh pardeshi

void cluster_means(int dim, int n, float *data, int k, int* cluster_assign, float* cluster_center, int* cluster_size, int cluster1, int cluster2);
void cluster_apportionment(int dim, int n, int k, float* data, float* cluster_center, int cluster1, int cluster2, int * cluster_assign);
float calc_distance(int dim, float* p1, float* p2, int p, int q);
void recompute_centroids(int dim, int n, int k, float* data, int* cluster_assign, int* cluster_size, int cluster1, int cluster2, float*);
void bisect_centroid_compute(int k, int** points, int* cluster_size, float* data, int dim, int one_cluster, int new_cluster, float*);
float* update_cluster_center(int dim, int n, float* data, float * cluster_center, int cluster_count);

float* calc_farthest_in_cluster(int dim, int cluster, float *data, int * cluster_assign, int n, float* cluster_center);
int bi_calc_ssd(int dim, int n, float* data, int cluster_count, float *cluster_center, int* cluster_assign);



float random()
{

	return (float)rand() / RAND_MAX;
}

float square_error(float* data, int* cluster_assign, float* cluster_center, int n, int k, int dim)
{
	int i, j;
	float distance = 0.0;
	for (i = 0; i<k; i++)
	{
		for (j = 0; j<n; j++)
		{
			if (cluster_assign[j] == i)
			{
				distance += sqrt(calc_distance(dim, data, cluster_center, j*dim, i*dim));
			}
		}
	}
	return distance;
}


void main()
{
	int i = 10000;
	int j = 32;
	int d = 3;
	int a, aa, bb, b;
	float *data = (float*)calloc(i*d, sizeof(float));
	int *bassign = (int*)calloc(i, sizeof(int));
	float *bcenter = (float*)calloc(j*d, sizeof(float));
	int *bsize = (int*)calloc(j, sizeof(int));
	float *bpoint = (float*)calloc(d, sizeof(float));
	float *bpoint1 = (float*)calloc(d, sizeof(float));
	int bcluster1 = 0;
	int bcluster2 = 1;
	int bcount = 1;
	int flag = 0;
	int val = 0;
	srand(1);
	for (a = 0; a<i; a++)
		for (b = 0; b<d; b++)
		{
			data[a*d + b] = random();
		}

	//choosing first center randomly
	memcpy(&bcenter[0], &data[0], d * sizeof(float));

	//kmeans clustering

	for (a = 0; a < i; a++)
		bassign[a] = 0;

	bsize[0] = i;

	bpoint = calc_farthest_in_cluster(d, bcluster1, data, bassign, i, bcenter);

	memcpy(&bcenter[d], &bpoint[0], d * sizeof(float));


	while (bcount != j)
	{
		cluster_means(d, i, data, j, bassign, bcenter, bsize, bcluster1, bcluster2);
		bcount++;
		
		bcluster1 = bi_calc_ssd(d, i, data, bcount, bcenter, bassign);

		for (aa = 0; aa < bcount; aa++)
		{
			if (bsize[aa] == 1)
			{
				flag = 1;
				bcount--;
				bpoint = update_cluster_center(d, i, data, bcenter, bcount);
				for (bb = 0; bb < d; bb++)
				{
					bcenter[aa*d + bb] = bpoint[bb];
				}

			}
		}
		bcluster2 = bcount;
		if (bcount == j)
			break;
		if (flag != 1)
		{
			while (1)
			{
				if (bassign[val] == bcluster1)
					break;

				val++;
			}

			memcpy(&bcenter[bcluster1*d], &data[val*d], d * sizeof(float));
			bpoint1 = calc_farthest_in_cluster(d, bcluster1, data, bassign, i, bcenter);
			memcpy(&bcenter[bcluster2*d], &bpoint1[0], d * sizeof(float));
		}

	}

	/*for (i = 0; i < k; i++)
	for (j = 0; j < dim; j++)
	{
	printf("Centers of cluster %d------Coordinate %d is %f \n", i, j, cluster_center[i*dim + j]);
	}*/

	for (a = 0; a < j; a++)
		printf("Size of cluster%d is %d\n", a, bsize[a]);

	float sse = square_error(data, bassign, bcenter, i, j, d);
	printf("SSE = %f\n", sse);


	getchar();

}






float* calc_farthest_in_cluster(int d, int bcluster, float *data, int * bassign, int i, float* bcenter)
{
	float *temp = (float*)calloc(d, sizeof(float));
	int a, b;
	float diff = 0.0;
	float max = 0.0;
	int bpoint;
	for (a = 0; a < i; a++)
	{
		if (bassign[a] == bcluster)
		{
			diff = calc_distance(d, data, bcenter, a*d, bcluster*d);
			if (max < diff)
			{
				max = diff;
				bpoint = i;
			}
		}
	}

	for (a = 0; a < d; a++)
	{
		temp[a] = data[bpoint*d + a];
	}

	return temp;

}


int bi_calc_ssd(int d, int i, float* data, int bcount, float *bcenter, int* bassign)
{

	int a, b;
	int c;
	float *val = (float*)calloc(bcount, sizeof(float));
	int max = 0;
	for (a = 0; a < bcount; a++)
	{
		for (b = 0; b < i; b++)
		{
			if (bassign[b] == a)
			{
				val[a] += calc_distance(d, data, bcenter, b*d, a*d);
			}
		}
		if (max < val[a])
		{
			max = val[a];
			c = a;
		}
	}

	return c;
}

float* update_cluster_center(int d, int i, float* data, float * bcenter, int bcount)
{
	int a, b;
	float diff;
	float max = 0.0;
	int x;
	float* difference = (float*)calloc(i, sizeof(float));
	float* temp = (float*)calloc(d, sizeof(float));
	for (a = 0; a < i; a++)
	{
		for (b = 0; b < bcount; b++)
		{
			diff = sqrt(calc_distance(d, data, bcenter, a*d, b*d));
			if (max < diff)
			{
				max = diff;
				difference[a] = diff;
			}
		}
		max = 0.0;
	}
	max = 0.0;
	for (a = 0; a < i; a++)
	{
		if (max < difference[a])
		{
			max = difference[a];
			x = a;
		}
	}
	for (a = 0; a < d; a++)
	{
		temp[a] = data[x*d + a];
	}

	return temp;
}


void cluster_means(int d, int i, float *data, int j, int* bassign, float* bcenter, int* bsize, int bcluster1, int bcluster2)
{

	float *center = (float*)calloc(2 * d, sizeof(float));
	int a;
	int aa, bb, cc;
	int rr;




	cluster_apportionment(d, i, j, data, bcenter, bcluster1, bcluster2, bassign);
	recompute_centroids(d, i, j, data, bassign, bsize, bcluster1, bcluster2, bcenter);

	for (aa = 0; aa < d; aa++)
	{
		bcenter[bcluster1*d + aa] = center[aa];
	}
	for (aa = 0; aa < d; aa++)
	{
		bcenter[bcluster2*d + aa] = center[d + aa];
	}
}


void cluster_apportionment(int d, int i, int j, float* data, float* bcenter, int bcluster1, int bcluster2, int * bassign)
{
	int a;
	int b;
	float min;
	int bcluster = 0;
	int *cluster_assigned = (int*)calloc(i, sizeof(int));
	float distance_from_center1 = 0.0;
	float distance_from_center2 = 0.0;

	for (a = 0; a < i; a++)
	{
		if (bassign[a] == bcluster1)
		{
			distance_from_center1 = sqrt(calc_distance(d, data, bcenter, a*d, bcluster1*d));
			distance_from_center2 = sqrt(calc_distance(d, data, bcenter, a*d, bcluster2*d));
			if (distance_from_center1 > distance_from_center2)
				bassign[a] = bcluster2;
			else
				bassign[a] = bcluster1;
		}
	}




}

float calc_distance(int d, float* x1, float* x2, int x, int y)
{

	int a;
	float value;
	float difference = 0;
	for (a = 0; a < d; a++)
	{

		difference += pow((x1[x + a] - x2[y + a]), 2);
	}
	return difference;
}

void recompute_centroids(int d, int i, int j, float* data, int* bassign, int* bsize, int bcluster1, int bcluster2, float* bcenter)
{
	int **bpoints = (int**)calloc(2, sizeof(int));	float *centers = (float*)calloc(2 * d, sizeof(float));
	int a;
	int b;
	int p1 = 0;
	int p2 = 0;

	for (a = 0; a < 2; a++)
	{
		bpoints[a] = (int*)calloc(i, sizeof(int));
	}


	for (b = 0; b < i; b++)
	{
		if (bassign[b] == bcluster1)
		{
			bpoints[0][p1++] = b;

		}
		else if (bassign[b] == bcluster2)
		{
			bpoints[1][p2++] = b;

		}

	}


	bsize[bcluster1] = p1;

	bsize[bcluster2] = p2;





	bisect_centroid_compute(j, bpoints, bsize, data, d, bcluster1, bcluster2, bcenter);



}

void bisect_centroid_compute(int j, int** bpoints, int* bsize, float* data, int d, int bcluster1, int bcluster2, float* bcenter)
{
	//float *centercentroids = (float*)calloc(2*dim, sizeof(float));
	int a;
	int b;
	int x = 0;
	int y = 0;
	float *coordinates = (float*)calloc(2 * d, sizeof(float));



	for (a = 0; a < d; a++)
	{
		while (x != bsize[bcluster1])
		{
			coordinates[a] += data[bpoints[0][x] * d + a];
			x++;
		}
		bcenter[bcluster1*d + a] = coordinates[a] / bsize[bcluster1];
		x = 0;
	}

	for (a = 0; a < d; a++)
		coordinates[a] = 0.0;

	for (a = 0; a < d; a++)
	{
		while (y != bsize[bcluster2])
		{
			coordinates[a] += data[bpoints[1][y] * d + a];
			y++;
		}
		bcenter[bcluster2*d + a] = coordinates[a] / bsize[bcluster2];
		y = 0;
	}


}