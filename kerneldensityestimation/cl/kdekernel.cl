

#include "samplers.cl" 
#include "image3d_write.cl" 

//__kernel void volumeKDEKernel(	global float* stencil, 
//								int3 stencil_dims, 
//								int3 stencil_half_dims,
//								int3 vol_dims, 							
//								global int* extremas,
//								int extrema_size,
//								image_3d_write_float32_t out_data	) 
//{
//	const int xy_ = vol_dims.x * vol_dims.y;
//	const float stencil_val = stencil[stencil_dims.x * stencil_dims.y * get_global_id(2) + stencil_dims.x * get_global_id(1) + get_global_id(0)];
//	for(int i = 0; i < extrema_size; ++i) {		
//		int3 extrema_idx = (int3) (extremas[i] % vol_dims.x, (int)(extremas[i] / vol_dims.x) % vol_dims.y, extremas[i] / xy_);
//		int3 global_idx = (int3) (	extrema_idx.x + get_global_id(0) - stencil_half_dims.x, 
//									extrema_idx.y + get_global_id(1) - stencil_half_dims.y, 
//									extrema_idx.z + get_global_id(2) - stencil_half_dims.z	);	
//		if(	global_idx.x > -1 && global_idx.x < vol_dims.x &&
//			global_idx.y > -1 && global_idx.y < vol_dims.y &&
//			global_idx.z > -1 && global_idx.z < vol_dims.z )
//		{
//			//out_data[0] += stencil_val;
//			out_data[xy_ * global_idx.z + vol_dims.x * global_idx.y + global_idx.x] += stencil_val;
//		}
//	}
//}

__kernel void volumeKDEKernel(	global float* stencil, 
								int3 stencil_dims, 
								int3 stencil_half_dims,
								int3 vol_dims, 							
								int3 extrema_idx,
								image_3d_write_float32_t out_data	) 
{
	const int xy_ = vol_dims.x * vol_dims.y;
	const float stencil_val = stencil[stencil_dims.x * stencil_dims.y * get_global_id(2) + stencil_dims.x * get_global_id(1) + get_global_id(0)];
	
	int3 global_idx = (int3) (	extrema_idx.x + get_global_id(0) - stencil_half_dims.x, 
								extrema_idx.y + get_global_id(1) - stencil_half_dims.y, 
								extrema_idx.z + get_global_id(2) - stencil_half_dims.z	);	
	if(	global_idx.x > -1 && global_idx.x < vol_dims.x &&
		global_idx.y > -1 && global_idx.y < vol_dims.y &&
		global_idx.z > -1 && global_idx.z < vol_dims.z )
	{
		//out_data[0] += stencil_val;
		out_data[xy_ * global_idx.z + vol_dims.x * global_idx.y + global_idx.x] += stencil_val;
	}
}
