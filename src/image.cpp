//Kutay Özbay
//270201017

#include "ceng391/image.hpp"

#include <cstdlib>
#include <cstdio>
#include<fstream>
#include<iostream>
#include<cstring>
#include <stdexcept>

#include <png.h>

#ifndef png_jmpbuf
#  define png_jmpbuf(png_ptr) ((png_ptr)->jmpbuf)
#endif

#ifndef png_infopp_NULL
#  define png_infopp_NULL (png_infopp)NULL
#endif

#ifndef int_p_NULL
# define int_p_NULL (int*)NULL
#endif

using namespace std;

namespace ceng391 {

Image::Image(int width, int height, int n_channels, int step)
{
        m_width = width;
        m_height = height;
        m_n_channels = n_channels;
        int row_len = m_width * m_n_channels;
        if (step < row_len)
                step = row_len;
        m_step = step;

        m_data = new uchar[m_step * m_height];
}

Image::~Image()
{
        delete [] m_data;
}

void Image::reallocate(int width, int height, int n_channels)
{
        if (width  != this->m_width ||
            height != this->m_height ||
            n_channels != this->m_n_channels) {
                delete [] m_data;
                int step = width * n_channels;
                m_step = step;
                m_data = new uchar[height * m_step];
                m_width = width;
                m_height = height;
                m_n_channels = n_channels;
        }
}

void Image::set_rect_rgba(int x_tl, int y_tl, int width, int height,
                          uchar red, uchar green, uchar blue, uchar alpha)
{
        if (x_tl < 0) {
                width += x_tl;
                x_tl = 0;
        }

        if (y_tl < 0) {
                height += y_tl;
                y_tl = 0;
        }

        size_t length = min(width, m_width - x_tl);
        int y_max = min(y_tl + height, m_height);
        if (m_n_channels == 4) {
                for (int y = y_tl; y < y_max; ++y) {
                        uchar *row_y = data(y) + x_tl;
                        for (int x = x_tl; x < x_tl + length; ++x) {
                                row_y[4*x]   = red;
                                row_y[4*x+1] = green;
                                row_y[4*x+2] = blue;
                                row_y[4*x+3] = alpha;
                        }
                }
        } else if (m_n_channels == 3) {
                for (int y = y_tl; y < y_max; ++y) {
                        uchar *row_y = data(y) + x_tl;
                        for (int x = x_tl; x < x_tl + length; ++x) {
                                row_y[3*x]   = red;
                                row_y[3*x+1] = green;
                                row_y[3*x+2] = blue;
                        }
                }
        } else if (m_n_channels == 1) {
                int value = 0.3f * red + 0.59f * green + 0.11f * blue;
                if (value > 255)
                        value = 255;
                for (int y = y_tl; y < y_max; ++y) {
                        uchar *row_y = data(y) + x_tl;
                        memset(row_y, value, length);
                }
        } else {
                cerr << "Unexpected number of channels ("
                     << m_n_channels << ") in call to set_rect_xxx()" << endl;
                exit(EXIT_FAILURE);
        }
}

void Image::xsave_pnm(const  std::string& filename) const
{
        if (m_n_channels != 1) {
                cerr << "xsave_pnm() works only with grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        const string magic_head = "P5";
        ofstream fout;
        string extended_name = filename + ".pgm";
        fout.open(extended_name.c_str(), ios::out | ios::binary);
        fout << magic_head << "\n";
        fout << m_width << " " << m_height << " 255\n";
        for (int y = 0; y < m_height; ++y) {
                const uchar *row_data = this->data(y);
                fout.write(reinterpret_cast<const char*>(row_data),
                           m_width*sizeof(uchar));
        }
        fout.close();
}

void Image::xsave_png(const std::string &filename) const
{
        // We open the output file with C style IO since we are using libpng
        // C-API
        FILE *fout = fopen(filename.c_str(), "w+b");
        if (!fout) {
                cerr << "[ERROR][CENG391::Image] Failed open file for writing: " << filename << endl;
                exit(EXIT_FAILURE);
        }

        png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                      0, 0, 0);
        if (!png_ptr) {
                cerr << "[ERROR][CENG391::Image] Failed to create PNG write structure!" << endl;
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
                cerr << "[ERROR][CENG391::Image] Failed to create PNG info structure!" << endl;
                png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        if (setjmp(png_jmpbuf(png_ptr))) {
                cerr << "[ERROR][CENG391::Image] Failed to create PNG jump buffer!" << endl;
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        int color_type;
        switch (this->m_n_channels) {
        case 1: color_type = PNG_COLOR_TYPE_GRAY; break;
        case 3: color_type = PNG_COLOR_TYPE_RGB; break;
        case 4: color_type = PNG_COLOR_TYPE_RGBA; break;
        default:
                cerr << "[ERROR][CENG391::Image] Unsupported image type for saving as PNG!" << endl;
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        png_init_io(png_ptr, fout);
        png_set_IHDR(png_ptr, info_ptr, this->m_width, this->m_height, 8,
                     color_type, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png_ptr, info_ptr);

        png_bytepp row_pointers = (png_bytepp)malloc(this->m_height * sizeof(png_bytep));
        if (!row_pointers) {
                cerr << "[ERROR][CENG391::Image]Error creating PNG row pointers" << endl;
                png_destroy_write_struct(&png_ptr, &info_ptr);
                fclose(fout);
                exit(EXIT_FAILURE);
        }

        for (png_int_32 k = 0; k < this->m_height; k++) {
                row_pointers[k] = (png_bytep)(this->data(k));
        }

        png_write_image(png_ptr, row_pointers);
        png_write_end(png_ptr, info_ptr);

        png_destroy_write_struct(&png_ptr, &info_ptr);
        free(row_pointers);
        fclose(fout);
}


void Image::xload_png(const std::string &filename, bool load_as_grayscale)
{
        // We open the output file with C style IO since we are using libpng
        // C-API
        FILE *fin = fopen(filename.c_str(), "r+b");
        if (!fin) {
                cerr << "[ERROR][CENG391::Image] Failed to open file for reading: " << filename << endl;
                exit(EXIT_FAILURE);
        }

        png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                                     NULL, NULL, NULL);
        if (!png_ptr) {
                cerr << "[ERROR][CENG391::Image] Could not create PNG read structure" << endl;
                fclose(fin);
                exit(EXIT_FAILURE);
        }

        png_infop info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) {
                cerr << "[ERROR][CENG391::Image] Could not create PNG info pointer" << endl;
                png_destroy_read_struct(&png_ptr, png_infopp_NULL,
                                        png_infopp_NULL);
                fclose(fin);
                exit(EXIT_FAILURE);
        }

        if (setjmp(png_jmpbuf(png_ptr))) {
                cerr << "[ERROR][CENG391::Image] Could not set jump point for reading PNG file" << endl;
                png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
                fclose(fin);
                exit(EXIT_FAILURE);
        }

        png_init_io(png_ptr, fin);
        png_read_info(png_ptr, info_ptr);

        png_uint_32 stream_width, stream_height;
        int bit_depth, color_type, interlace_type;
        png_get_IHDR(png_ptr, info_ptr, &stream_width, &stream_height, &bit_depth, &color_type,
                     &interlace_type, int_p_NULL, int_p_NULL);

        png_set_strip_16(png_ptr);
        if (color_type == PNG_COLOR_TYPE_GA) {
                png_set_strip_alpha(png_ptr); /*(not recommended). */
        }

        png_set_packing(png_ptr);
        if (color_type == PNG_COLOR_TYPE_PALETTE) {
                png_set_palette_to_rgb(png_ptr);
        }

        png_set_expand(png_ptr);

        // Depending on the type of image in the file and the load_as_grayscale
        // flag, we determine the desired number of channels of the output
        // image.
        int desired_n_channels = 4;
        if (load_as_grayscale) {
                desired_n_channels = 1;
                png_set_rgb_to_gray_fixed(png_ptr, 1, 30000, 59000);
                png_set_strip_alpha(png_ptr); /*(not recommended). */
        } else {
                if (color_type == PNG_COLOR_TYPE_GRAY ||
                    color_type == PNG_COLOR_TYPE_GA) {
                        desired_n_channels = 1;
                }

                if(color_type == PNG_COLOR_TYPE_RGB ||
                   color_type == PNG_COLOR_TYPE_PALETTE) {
                        png_set_add_alpha(png_ptr, 255, PNG_FILLER_AFTER);
                }
        }

        // If the current image dimensions do not match the image to be loaded,
        // then reallocate with the desired dimensions.
        reallocate(stream_width, stream_height, desired_n_channels);

        png_bytepp row_pointers = (png_bytepp)malloc(this->m_height * sizeof(png_bytep));
        if (!row_pointers) {
                cerr << "[ERROR][CENG391::Image]Error creating PNG row pointers" << endl;
                png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
                fclose(fin);
                exit(EXIT_FAILURE);
        }
        for (int k = 0; k < this->m_height; k++) {
                row_pointers[k] = (png_bytep)(this->data(k));
        }

        png_read_image(png_ptr, row_pointers);
        png_read_end(png_ptr, info_ptr);

        png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);

        free(row_pointers);

        fclose(fin);
}

void copy_row_to_buffer(uchar* buffer, const uchar* row,
                        int n, int border_size)
{
        memcpy(reinterpret_cast<void *>(buffer + border_size),
               reinterpret_cast<const void *>(row),
               n*sizeof(*buffer));

        for (int i = 0; i < border_size; ++i) {
                buffer[i] = buffer[border_size];
                buffer[n-i-1] = buffer[n+border_size-1];
        }
}

void copy_column_to_buffer(uchar* buffer, const uchar* column,
                           int n, int step, int border_size)
{
        for (int i = 0; i < n; ++i) {
                buffer[border_size + i] = column[i * step];
        }

        for (int i = 0; i < border_size; ++i) {
                buffer[i] = buffer[border_size];
                buffer[n-i-1] = buffer[n+border_size-1];
        }
}

void box_filter_buffer(int n, uchar* buffer, int filter_size)
{
        for(int i = 0; i < n; ++i) {
                int sum = 0;
                for (int j = 0; j < filter_size; ++j)
                        sum += buffer[i + j];
                sum /= filter_size;
                buffer[i] = sum;
        }
}

void Image::box_filter_x(int filter_size)
{
        if (m_n_channels != 1) {
                cerr << "Can only box filter grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int border_size = filter_size / 2;
        filter_size = 2*border_size + 1;

        int lbuffer = 2*border_size + m_width;
        uchar *buffer = new uchar[lbuffer];

        for (int y = 0; y < m_height; ++y) {
                copy_row_to_buffer(buffer, data(y), m_width, border_size);
                box_filter_buffer(m_width, buffer, filter_size);
                memcpy(reinterpret_cast<void *>(data(y)),
                       reinterpret_cast<const void *>(buffer),
                       m_width*sizeof(*buffer));
        }

        delete [] buffer;
}

void Image::box_filter_y(int filter_size)
{
        if (m_n_channels != 1) {
                cerr << "Can only box filter grayscale images" << endl;
                exit(EXIT_FAILURE);
        }
        int border_size = filter_size / 2;
        filter_size = 2*border_size + 1;

        int lbuffer = 2*border_size + m_height;
        uchar *buffer = new uchar[lbuffer];

        for (int x = 0; x < m_width; ++x) {
                copy_column_to_buffer(buffer, m_data + x, m_height,
                                      m_step, border_size);
                box_filter_buffer(m_height, buffer, filter_size);
                for (int y = 0; y < m_height; ++y)
                        m_data[x + y*m_step] = buffer[y];
        }

        delete [] buffer;
}

void Image::box_filter(int filter_size)
{
        box_filter_x(filter_size);
        box_filter_y(filter_size);
}



int accessPixel(unsigned char * arr, int col, int row, int k, int width, int height) 
{
    int sum = 0;
    int sumKernel = 0;

    for (int j = -1; j <= 1; j++) 
    {
        for (int i = -1; i <= 1; i++) 
        {
            if ((row + j) >= 0 && (row + j) < height && (col + i) >= 0 && (col + i) < width) 
            {
                int color = arr[(row + j) * 3 * width + (col + i) * 3 + k];
                sum += color * kernel[i + 1][j + 1];
                sumKernel += kernel[i + 1][j + 1];
            }
        }
    }

    return sum / sumKernel;
}
}
Image Image::smooth_x(float sigma_x, int l){
        
    double r, s = 2.0 * sigma_x * sigma_x;
 
    // sum is for normalization
    double sum = 0.0;

    filter_length = (2*l) + 1;
 
    // generating 5x5 kernel
    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[x + 2][y + 2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[x + 2][y + 2];
        }
    }
 
    // normalising the Kernel
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            GKernel[i][j] /= sum;


    for (int row = 0; row < img.h(); row++) 
    {
        for (int col = 0; col < img.w(); col++) 
        {
            for (int k = 0; k < 3; k++) 
            {
                result[3 * row * width + 3 * col + k] = accessPixel(img, col, row, k, img.w(), img.h());
            }
        }
    }


}

void Image::smooth_y(float sigma_y){
        
    double r, s = 2.0 * sigma_y * sigma_y;
 
    // sum is for normalization
    double sum = 0.0;

    filter_length = (2*l) + 1;
 
    // generating 5x5 kernel
    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[x + 2][y + 2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[x + 2][y + 2];
        }
    }
 
    // normalising the Kernel
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            GKernel[i][j] /= sum;

    for (int row = 0; row < img.h(); row++) 
    {
        for (int col = 0; col < img.w(); col++) 
        {
            for (int k = 0; k < 3; k++) 
            {
                result[3 * row * width + 3 * col + k] = accessPixel(img, col, row, k, img.w(), img.h());
            }
        }
    }
}

void Image::smooth(float sigma_x, float sigma_y){
        
    double r, s = 2.0 * sigma_x * sigma_y;
 
    // sum is for normalization
    double sum = 0.0;

    filter_length = (2*l) + 1;
 
    // generating 5x5 kernel
    for (int x = -2; x <= 2; x++) {
        for (int y = -2; y <= 2; y++) {
            r = sqrt(x * x + y * y);
            GKernel[x + 2][y + 2] = (exp(-(r * r) / s)) / (M_PI * s);
            sum += GKernel[x + 2][y + 2];
        }
    }
 
    // normalising the Kernel
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            GKernel[i][j] /= sum;

    for (int row = 0; row < img.h(); row++) 
    {
        for (int col = 0; col < img.w(); col++) 
        {
            for (int k = 0; k < 3; k++) 
            {
                result[3 * row * width + 3 * col + k] = accessPixel(img, col, row, k, img.w(), img.h());
            }
        }
    }
}

void Image::deriv_x(){
        int a[3][3] = {  
        {-1, 0, 1} ,   /*  initializers for row indexed by 0 */
        {-2, 0, 2} ,   /*  initializers for row indexed by 1 */
        {-1, 0, 1}   /*  initializers for row indexed by 2 */
        };
        // Handle border issues
        Mat1b _src;
        copyMakeBorder(src, _src, radius, radius, radius, radius, BORDER_REFLECT101);

        // Create output matrix
        dst.create(src.rows, src.cols);

        // Convolution loop

        // Iterate on image 
        for (int r = radius; r < _src.rows - radius; ++r)
        {
                for (int c = radius; c < _src.cols - radius; ++c)
                {
                        short s = 0;

                        // Iterate on kernel
                        for (int i = -radius; i <= radius; ++i)
                        {
                                for (int j = -radius; j <= radius; ++j)
                                {
                                        s += _src(r + i, c + j) * kernel(i + radius, j + radius);
                                }
                        }
                        dst(r - radius, c - radius) = s;
                }       
        }
}

void Image::deriv_y(){
        int a[3][3] = {  
        {-1, -2, -1} ,   /*  initializers for row indexed by 0 */
        {0, 0, 0} ,   /*  initializers for row indexed by 1 */
        {1, 2, 1}   /*  initializers for row indexed by 2 */
        };
                // Handle border issues
        Mat1b _src;
        copyMakeBorder(src, _src, radius, radius, radius, radius, BORDER_REFLECT101);

        // Create output matrix
        dst.create(src.rows, src.cols);

        // Convolution loop

        // Iterate on image 
        for (int r = radius; r < _src.rows - radius; ++r)
        {
                for (int c = radius; c < _src.cols - radius; ++c)
                {
                        short s = 0;

                        // Iterate on kernel
                        for (int i = -radius; i <= radius; ++i)
                        {
                                for (int j = -radius; j <= radius; ++j)
                                {
                                        s += _src(r + i, c + j) * kernel(i + radius, j + radius);
                                }
                        }
                        dst(r - radius, c - radius) = s;
                }       
        }
}
}

Image Image::rotate(Image img, float theta){


 int row1=0,col1=0,row2=0,col2=0;

 double aff[3][3];
 aff[0][0]=cos(theta);
 aff[0][1]=-sin(theta);
 aff[0][2]=0;

 aff[1][0]=sin(theta);
 aff[1][1]=cos(theta);
 aff[1][2]=0;

 aff[2][0]=0;
 aff[2][1]=0;
 aff[2][2]=1;



 int mx=(img.cols/2);
 int my=(img.rows/2);
 for(int i=0;i<img.rows;i++) // y
    {
      for (int j=0;j<img.cols;j++) // x
      {

           int x=(j-mx) *aff[0][0]+(i-my)*aff[0][1]+aff[0][2];

           int y=(j-mx) *aff[1][0]+(i-my)*aff[1][1]+aff[1][2];

           row1=min(row1,y); row2=max(row2,y);
           col1=min(col1,x); col2=max(col2,x);

      }
    }

 int dg=img.rows;
 int dgx=img.cols;

 Image nimg(dg,dgx, CV_8U, Scalar(0,0,0));


 for(int i=0;i<img.rows;i++) // y
 {
      for (int j=0;j<img.cols;j++) // x
      {
        nimg.at<uchar>(i+(nimg.rows-img.rows)/2,j+(nimg.cols-img.cols)/2) =img.at<uchar>(i,j);

      }
 }

 int midx=(nimg.cols/2);
 int midy=(nimg.rows/2);
 Mat img(dg,dgx, CV_8U, Scalar(0,0,0));

 double mdx=(img.cols/2);
 double mdy=(img.rows/2);


 for(int i=0;i<img.rows;i++) // y
 {
      for (int j=0;j<img.cols;j++) // x
      {

        double x=(j-mdx)*aff[0][0]+(i-mdy)*aff[1][0]-(aff[0][2]*aff[0][0]+aff[1][2]*aff[1][0]);

        double y=(j-mdx)*aff[0][1]+(i-mdy)*aff[1][1]+-(aff[0][2]*aff[0][1]+aff[1][2]*aff[1][1]);

        if ( (int)y+midy>=0 && (int)y+midy<nimg.cols && (int)x+midx>=0 && (int)x+midx<nimg.rows  )
              {
                if (bilinear==1){
                int x1=(int)xx; int x2=x1+1;
                int y1=(int) yy; int y2=y1+1;


                vector< pair<int,int> > nb;

                nb.push_back(make_pair(x1,y1));
                nb.push_back(make_pair(x2,y1));
                nb.push_back(make_pair(x2,y2));
                nb.push_back(make_pair(x1,y2));

                double sum=0;
                int f=0;
                for (int i=0;i<nb.size();i++)
	        {

	        int x=nb[i].first,y=nb[i].second;
	        if ( (x<0 || x>=img.cols) || (y<0 || y>=img.rows))
		{}
                }
                else
                {  
                        double v=0;
                        v=  (1-abs(xx-x) )* (1- abs(yy-y))  ;
                        sum+=v*im(y,x) ;
                }
                        
                im(i,j)= sum;
	        }
                
        else
              im(i,j)=imo((int)y+midx,(int)x+midy);
        }


      }
 }


 return img;
 // end affine

}
