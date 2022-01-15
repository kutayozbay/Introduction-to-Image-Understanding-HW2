//Kutay Özbay
//270201017


#include <cstdlib>
#include <iostream>

#include "ceng391/image.hpp"

using namespace std;
using ceng391::Image;

int main(int argc, char** argv)
{
        Image *img = nullptr;
        if (argc == 1) {
                img = Image::new_grayscale(320, 240);
                cout << "Created image of size (" << img->w()
                     << "x" << img->h() << ")!" << endl;
                img->set(255);
        } else {
                img = Image::from_png(argv[1], false);
        }
        img->set_rect_rgba(100, -50, 1000, 150, 0, 255, 0, 128);
        img->xsave_pnm("/tmp/test_img");
        img->xsave_png("/tmp/test_img.png");
        delete img;
        return EXIT_SUCCESS;
}

/// Local Variables:
/// mode: c++
/// compile-command: "make -C ../build image-test"
/// End:
