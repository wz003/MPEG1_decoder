//
//  main.cpp
//  MPEG1
//
//  Created by 至 on 2017/6/19.
//  Copyright © 2017年 至. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>
using namespace std;

int frame_num=0;

const char hex_chars[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };

const int scan[8][8]={0,1,5,6,14,15,27,28,2,4,7,13,16,26,29,42,3,8,12,17,25,30,41,43,9,11,18,24,31,40,44,53,10,19,23,32,39,45,52,54,20,22,33,38,46,51,55,60,21,34,37,47,50,56,59,61,35,36,48,49,57,58,62,63};
const int intra_quan[8][8]={8,16,19,22,26,27,29,34,
                            16,16,22,24,27,29,34,37,
                            19,22,26,27,29,34,34,38,
                            22,22,26,27,29,34,37,40,
                            22,26,27,29,32,35,40,48,
                            26,27,29,32,35,40,48,58,
                            26,27,29,34,38,46,56,69,
                            27,29,35,38,46,56,69,83};

//macro_address_increment
map<string,int>macro_address;
//macroblock_type
map<string,string>macroblock_type_I;
map<string,string>macroblock_type_P;
map<string,string>macroblock_type_B;
map<string,string>macroblock_type_D;
//motion vector
map<string,int>macroblock_pattern;
map<string,int>motion_vector;
//
//DCT coe
map<string,int>dct_luminance;
map<string,int>dct_chrominance;
map<string,int>dct_coeff_run;
map<string,int>dct_coeff_level;

///

int tsn;// temperal sequence number
int ft; // frame type

int quan;
int pre_mb_addr;
int past_intra_addr;
int dct_dc_y_past;
int dct_dc_cb_past;
int dct_dc_cr_past;


int forward_f;
int forward_r_size;
int backward_f;
int backward_r_size;
int motion_horizontal_forward_r;
int motion_vertical_forward_r;
int motion_horizontal_backward_r;
int motion_vertical_backward_r;
int motion_horizontal_forward_code;
int motion_vertical_forward_code;
int motion_horizontal_backward_code;
int motion_vertical_backward_code;

int recon_right_for_prev;
int recon_down_for_prev;
int recon_right_for;
int recon_down_for;

int recon_right_back_prev;
int recon_down_back_prev;

int full_pel_forward_vector;
int full_pel_backward_vector;

string past_block_type;
int past_for_r;
int past_for_d;
int past_back_r;
int past_back_d;

int all_recon[300][6][8][8];

int ref_old=0;
int ref_new=0;

int help=0;

int frame_count=0;

void IDCT2(double *F, double *f)
{
    const double pi = acos(-1);
    const double c30 = (2*cos(pi/8));
    const double c45 = sqrt(2);
    const double c60 = (2*sin(pi/8));
    const double Q = c30 - c60;
    const double R = c30 + c60;
    double tmp1;
    double m[8];
    
    m[0] = F[0];
    m[1] = F[4];
    m[2] = F[2] - F[6];
    m[3] = F[2] + F[6];
    m[4] = F[5] - F[3];
    m[5] = F[1] + F[7] - F[3] - F[5];
    m[6] = F[1] - F[7];
    m[7] = F[1] + F[7] + F[3] + F[5];
    
    tmp1 = c60 * (m[4] + m[6]);
    F[0] = m[0];
    F[1] = m[1];
    F[2] = m[2] * c45;
    F[3] = m[3];
    F[4] = -Q * m[4] - tmp1;
    F[5] = m[5] * c45;
    F[6] = R * m[6] - tmp1;
    F[7] = m[7];
    
    m[0] = F[6] - F[7] - F[5];
    m[1] = F[0] - F[1];
    m[2] = F[2] - F[3];
    m[3] = F[0] + F[1];
    m[4] = F[6] - F[7];
    m[5] = F[4];
    m[6] = F[3];
    m[7] = F[7];
    
    F[0] = m[7];
    F[1] = m[0];
    F[2] = m[4];
    F[3] = m[1] + m[2];
    F[4] = m[3] + m[6];
    F[5] = m[1] - m[2];
    F[6] = m[3] - m[6];
    F[7] = m[5] - m[0];
    
    f[0*8] = F[4] + F[0];
    f[1*8] = F[3] + F[2];
    f[2*8] = F[5] - F[1];
    f[3*8] = F[6] - F[7];
    f[4*8] = F[6] + F[7];
    f[5*8] = F[5] + F[1];
    f[6*8] = F[3] - F[2];
    f[7*8] = F[4] - F[0];
}

void IDCT(int a[8][8])
{
    const double pi = acos(-1);
    const double X[8] =
    {
        sqrt(2)*cos(0 *pi/16)/4,
        cos(1*pi/16)/2,
        cos(2*pi/16)/2,
        cos(3*pi/16)/2,
        cos(4*pi/16)/2,
        cos(5*pi/16)/2,
        cos(6*pi/16)/2,
        cos(7*pi/16)/2
    };
    
    int i, j;
    double x[64], y[64];
    
    for(i = 0; i < 8; i++)
        for(j = 0; j < 8; j++)
            x[i*8+j] = a[i][j] * X[i] * X[j];
    
    for(i = 0; i < 8; i++)
        IDCT2(x+8*i, y+i);
    
    for(i = 0; i < 8; i++)
        IDCT2(y+8*i, x+i);
    
    for(i = 0; i < 8; i++)
        for(j = 0; j < 8; j++)
            a[i][j] = (int)(x[i*8+j] + 0.5);
}

int gop[20][240][320][3];

void group_write()
{
    
    int width=320;
    int height=240;
    
    
    
    typedef struct                       /**** BMP file header structure ****/
    {
        unsigned int   bfSize;           /* Size of file */
        unsigned short bfReserved1;      /* Reserved */
        unsigned short bfReserved2;      /* ... */
        unsigned int   bfOffBits;        /* Offset to bitmap data */
    } BITMAPFILEHEADER;
    
    typedef struct                       /**** BMP file info structure ****/
    {
        unsigned int   biSize;           /* Size of info header */
        int            biWidth;          /* Width of image */
        int            biHeight;         /* Height of image */
        unsigned short biPlanes;         /* Number of color planes */
        unsigned short biBitCount;       /* Number of bits per pixel */
        unsigned int   biCompression;    /* Type of compression to use */
        unsigned int   biSizeImage;      /* Size of image data */
        int            biXPelsPerMeter;  /* X pixels per meter */
        int            biYPelsPerMeter;  /* Y pixels per meter */
        unsigned int   biClrUsed;        /* Number of colors used */
        unsigned int   biClrImportant;   /* Number of important colors */
    } BITMAPINFOHEADER;
    
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    
    /* Magic number for file. It does not fit in the header structure due to alignment requirements, so put it outside */
    unsigned short bfType=0x4d42;
    bfh.bfReserved1 = 0;
    bfh.bfReserved2 = 0;
    bfh.bfSize = 2+sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)+width*height;
    bfh.bfOffBits = 0x36;
    
    bih.biSize = sizeof(BITMAPINFOHEADER);
    bih.biWidth = width;
    bih.biHeight = height;
    bih.biPlanes = 1;
    bih.biBitCount = 24;
    bih.biCompression = 0;
    bih.biSizeImage = 0;
    bih.biXPelsPerMeter = 5000;
    bih.biYPelsPerMeter = 5000;
    bih.biClrUsed = 0;
    bih.biClrImportant = 0;
    
    
    
    for (int k=0;k<frame_num-frame_count;k++){
        
        int pic[240][320][3]={0};
        int temp_RGB[3];
        for(int i=0;i<240;i++){
            for (int j=0;j<320;j++) {
                int y = gop[k][i][j][0];
                int cb = gop[k][i][j][1]-128;
                int cr = gop[k][i][j][2]-128;
                temp_RGB[0] = y+ 1.402*cr ;//R
                temp_RGB[1] = y- 0.34414*cb - 0.71414*cr ;//G
                temp_RGB[2] = y+ 1.772*cb; //B
                for (int l = 0;l<3;l++) {
                    if (temp_RGB[l]>255) {
                        temp_RGB[l]=255;
                    }else if (temp_RGB[l]<0){
                        temp_RGB[l]=0;
                    }
                    pic[i][j][l] = temp_RGB[l];
                }
            }
        }
        
        char filename[20];
        sprintf(filename,"frame%d.bmp", frame_count+k+1);
        
        FILE *file = fopen(filename, "wb");
        if (!file)
        {
            printf("Could not write file\n");
            return;
        }
        
        /*Write headers*/
        fwrite(&bfType,1,sizeof(bfType),file);
        fwrite(&bfh, 1, sizeof(bfh), file);
        fwrite(&bih, 1, sizeof(bih), file);
        
        /*Write bitmap*/
        for (int y = bih.biHeight-1; y>=0; y--) /*Scanline loop backwards*/
        {
            for (int x = 0; x < bih.biWidth; x++) /*Column loop forwards*/
            {
                /*compute some pixel values*/
                unsigned char r = pic[y][x][0];
                unsigned char g = pic[y][x][1];
                unsigned char b = pic[y][x][2];
                fwrite(&b, 1, 1, file);
                fwrite(&g, 1, 1, file);
                fwrite(&r, 1, 1, file);
            }
        }
        fclose(file);
        
    }
    frame_count=frame_num;
}

void imgwrite()
{
    
    int w_b =20;
    int width=320;
    int height=240;
    
    //
    //    //reconstruct Macroblock
    int MB[300][16][16][3];
    for (int j=0;j<300;j++){
        for(int m = 0;m<8;m++){
            for (int n =0;n<8;n++) {
                MB[j][m][n][0] = all_recon[j][0][m][n];
                MB[j][m][n+8][0] = all_recon[j][1][m][n];
                MB[j][m+8][n][0] = all_recon[j][2][m][n];
                MB[j][m+8][n+8][0] = all_recon[j][3][m][n];
                
                MB[j][2*m][2*n][1] = all_recon[j][4][m][n];
                MB[j][2*m][2*n+1][1] = all_recon[j][4][m][n];
                MB[j][2*m+1][2*n][1] = all_recon[j][4][m][n];
                MB[j][2*m+1][2*n+1][1] = all_recon[j][4][m][n];
                
                MB[j][2*m][2*n][2] = all_recon[j][5][m][n];
                MB[j][2*m][2*n+1][2] = all_recon[j][5][m][n];
                MB[j][2*m+1][2*n][2] = all_recon[j][5][m][n];
                MB[j][2*m+1][2*n+1][2] = all_recon[j][5][m][n];
            }
        }
    }
    
    for (int j=0;j<300;j++) {
        int hor = j%w_b;
        int ver = j/w_b;
        //printf("now:%d,hor:%d,ver:%d\n",j,hor,ver);
        for (int k=0;k<16;k++) {
            for (int l=0;l<16;l++) {
                for (int a=0;a<3;a++){
                    if (MB[j][k][l][a]>255) MB[j][k][l][a]=255;
                    if (MB[j][k][l][a]<0) MB[j][k][l][a]=0;
                }
                gop[tsn][ver*16+k][hor*16+l][0] = MB[j][k][l][0];
                gop[tsn][ver*16+k][hor*16+l][1] = MB[j][k][l][1];
                gop[tsn][ver*16+k][hor*16+l][2] = MB[j][k][l][2];
            }
        }
    }
    //printf("pic location:%d\n",tsn);
    int pic[240][320][3]={0};
    int temp_RGB[3];
    for(int i=0;i<240;i++){
        for (int j=0;j<320;j++) {
            int y = gop[tsn][i][j][0];
            int cb = gop[tsn][i][j][1]-128;
            int cr = gop[tsn][i][j][2]-128;
            temp_RGB[0] = y+ 1.402*cr ;//R
            temp_RGB[1] = y- 0.34414*cb - 0.71414*cr ;//G
            temp_RGB[2] = y+ 1.772*cb; //B
            for (int l = 0;l<3;l++) {
                if (temp_RGB[l]>255) {
                    temp_RGB[l]=255;
                }else if (temp_RGB[l]<0){
                    temp_RGB[l]=0;
                }
                pic[i][j][l] = temp_RGB[l];
            }
        }
    }
    

    
    typedef struct                       /**** BMP file header structure ****/
    {
        unsigned int   bfSize;           /* Size of file */
        unsigned short bfReserved1;      /* Reserved */
        unsigned short bfReserved2;      /* ... */
        unsigned int   bfOffBits;        /* Offset to bitmap data */
    } BITMAPFILEHEADER;
    
    typedef struct                       /**** BMP file info structure ****/
    {
        unsigned int   biSize;           /* Size of info header */
        int            biWidth;          /* Width of image */
        int            biHeight;         /* Height of image */
        unsigned short biPlanes;         /* Number of color planes */
        unsigned short biBitCount;       /* Number of bits per pixel */
        unsigned int   biCompression;    /* Type of compression to use */
        unsigned int   biSizeImage;      /* Size of image data */
        int            biXPelsPerMeter;  /* X pixels per meter */
        int            biYPelsPerMeter;  /* Y pixels per meter */
        unsigned int   biClrUsed;        /* Number of colors used */
        unsigned int   biClrImportant;   /* Number of important colors */
    } BITMAPINFOHEADER;
    
    BITMAPFILEHEADER bfh;
    BITMAPINFOHEADER bih;
    
    /* Magic number for file. It does not fit in the header structure due to alignment requirements, so put it outside */
    unsigned short bfType=0x4d42;
    bfh.bfReserved1 = 0;
    bfh.bfReserved2 = 0;
    bfh.bfSize = 2+sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER)+width*height;
    bfh.bfOffBits = 0x36;
    
    bih.biSize = sizeof(BITMAPINFOHEADER);
    bih.biWidth = width;
    bih.biHeight = height;
    bih.biPlanes = 1;
    bih.biBitCount = 24;
    bih.biCompression = 0;
    bih.biSizeImage = 0;
    bih.biXPelsPerMeter = 5000;
    bih.biYPelsPerMeter = 5000;
    bih.biClrUsed = 0;
    bih.biClrImportant = 0;
    
    FILE *file = fopen("output.bmp", "wb");
    if (!file)
    {
        printf("Could not write file\n");
        return;
    }
    
    /*Write headers*/
    fwrite(&bfType,1,sizeof(bfType),file);
    fwrite(&bfh, 1, sizeof(bfh), file);
    fwrite(&bih, 1, sizeof(bih), file);
    
    /*Write bitmap*/
    for (int y = bih.biHeight-1; y>=0; y--) /*Scanline loop backwards*/
    {
        for (int x = 0; x < bih.biWidth; x++) /*Column loop forwards*/
        {
            /*compute some pixel values*/
            unsigned char r = pic[y][x][0];
            unsigned char g = pic[y][x][1];
            unsigned char b = pic[y][x][2];
            fwrite(&b, 1, 1, file);
            fwrite(&g, 1, 1, file);
            fwrite(&r, 1, 1, file);
        }
    }
    fclose(file);
    

}

void build_macro_address()
{
    macro_address["1"]=1;
    macro_address["011"]=2;
    macro_address["010"]=3;
    macro_address["0011"]=4;
    macro_address["0010"]=5;
    macro_address["00011"]=6;
    macro_address["00010"]=7;
    macro_address["0000111"]=8;
    macro_address["0000110"]=9;
    macro_address["00001011"]=10;
    macro_address["00001010"]=11;
    macro_address["00001001"]=12;
    macro_address["00001000"]=13;
    macro_address["00000111"]=14;
    macro_address["00000110"]=15;
    macro_address["0000010111"]=16;
    macro_address["0000010110"]=17;
    macro_address["0000010101"]=18;
    macro_address["0000010100"]=19;
    macro_address["0000010011"]=20;
    macro_address["0000010010"]=21;
    macro_address["00000100011"]=22;
    macro_address["00000100010"]=23;
    macro_address["00000100001"]=24;
    macro_address["00000100000"]=25;
    macro_address["00000011111"]=26;
    macro_address["00000011110"]=27;
    macro_address["00000011101"]=28;
    macro_address["00000011100"]=29;
    macro_address["00000011011"]=30;
    macro_address["00000011010"]=31;
    macro_address["00000011001"]=32;
    macro_address["00000011000"]=33;
    macro_address["00000001111"]=34;//macroblock_stuffing
    macro_address["00000001000"]=35;//macroblock_escape
}

void build_macroblock_type()
{
    //1 : macroblock_ quant
    //2 : macroblock_ motion_ forward
    //3 : macroblock motion_ backward
    //4 : macroblock_ pattern
    //5 : macroblock_ intra
    
    
    //I
    macroblock_type_I["1"]="00001";
    macroblock_type_I["01"]="10001";
    //P
    macroblock_type_P["1"]="01010";
    macroblock_type_P["01"]="00010";
    macroblock_type_P["001"]="01000";
    macroblock_type_P["00011"]="00001";
    macroblock_type_P["00010"]="11010";
    macroblock_type_P["00001"]="10010";
    macroblock_type_P["000001"]="10001";
    //B
    macroblock_type_B["10"]="01100";
    macroblock_type_B["11"]="01110";
    macroblock_type_B["010"]="00100";
    macroblock_type_B["011"]="00110";
    macroblock_type_B["0010"]="01000";
    macroblock_type_B["0011"]="01010";
    macroblock_type_B["00011"]="00001";
    macroblock_type_B["00010"]="11110";
    macroblock_type_B["000011"]="11010";
    macroblock_type_B["000010"]="10110";
    macroblock_type_B["000001"]="10001";
    //D
    macroblock_type_D["1"]="00001";
    
}
void build_macroblock_pattern()
{
    macroblock_pattern["111"]=60;
    macroblock_pattern["1101"]=4;
    macroblock_pattern["1100"]=8;
    macroblock_pattern["1011"]=16;
    macroblock_pattern["1010"]=32;
    
    macroblock_pattern["10011"]=12;
    macroblock_pattern["10010"]=48;
    macroblock_pattern["10001"]=20;
    macroblock_pattern["10000"]=40;
    macroblock_pattern["01111"]=28;
    
    macroblock_pattern["01110"]=44;
    macroblock_pattern["01101"]=52;
    macroblock_pattern["01100"]=56;
    macroblock_pattern["01011"]=1;
    macroblock_pattern["01010"]=61;
    
    macroblock_pattern["01001"]=2;
    macroblock_pattern["01000"]=62;
    macroblock_pattern["001111"]=24;
    macroblock_pattern["001110"]=36;
    macroblock_pattern["001101"]=3;
    
    macroblock_pattern["001100"]=63;
    macroblock_pattern["0010111"]=5;
    macroblock_pattern["0010110"]=9;
    macroblock_pattern["0010101"]=17;
    macroblock_pattern["0010100"]=33;
    
    macroblock_pattern["0010011"]=6;
    macroblock_pattern["0010010"]=10;
    macroblock_pattern["0010001"]=18;
    macroblock_pattern["0010000"]=34;
    macroblock_pattern["00011111"]=7;
    
    macroblock_pattern["00011110"]=11;
    macroblock_pattern["00011101"]=19;
    
    macroblock_pattern["00011100"]=35;
    macroblock_pattern["00011011"]=13;
    macroblock_pattern["00011010"]=49;
    macroblock_pattern["00011001"]=21;
    macroblock_pattern["00011000"]=41;
    
    macroblock_pattern["00010111"]=14;
    macroblock_pattern["00010110"]=50;
    macroblock_pattern["00010101"]=22;
    macroblock_pattern["00010100"]=42;
    macroblock_pattern["00010011"]=15;
    
    macroblock_pattern["00010010"]=51;
    macroblock_pattern["00010001"]=23;
    macroblock_pattern["00010000"]=43;
    macroblock_pattern["00001111"]=25;
    macroblock_pattern["00001110"]=37;
    
    macroblock_pattern["00001101"]=26;
    macroblock_pattern["00001100"]=38;
    macroblock_pattern["00001011"]=29;
    macroblock_pattern["00001010"]=45;
    macroblock_pattern["00001001"]=53;
    
    macroblock_pattern["00001000"]=57;
    macroblock_pattern["00000111"]=30;
    macroblock_pattern["00000110"]=46;
    macroblock_pattern["00000101"]=54;
    macroblock_pattern["00000100"]=58;
    
    macroblock_pattern["000000111"]=31;
    macroblock_pattern["000000110"]=47;
    macroblock_pattern["000000101"]=55;
    macroblock_pattern["000000100"]=59;
    macroblock_pattern["000000011"]=27;
    
    macroblock_pattern["000000010"]=39;
}
void build_motion_vector()
{
    motion_vector["00000011001"]=-16;
    motion_vector["00000011011"]=-15;
    motion_vector["00000011101"]=-14;
    motion_vector["00000011111"]=-13;
    motion_vector["00000100001"]=-12;
    motion_vector["00000100011"]=-11;
    motion_vector["0000010011"]=-10;
    motion_vector["0000010101"]=-9;
    motion_vector["0000010111"]=-8;
    motion_vector["00000111"]=-7;
    motion_vector["00001001"]=-6;
    motion_vector["00001011"]=-5;
    motion_vector["0000111"]=-4;
    motion_vector["00011"]=-3;
    motion_vector["0011"]=-2;
    motion_vector["011"]=-1;
    motion_vector["1"]=0;
    motion_vector["010"]=1;
    motion_vector["0010"]=2;
    motion_vector["00010"]=3;
    motion_vector["0000110"]=4;
    motion_vector["00001010"]=5;
    motion_vector["00001000"]=6;
    motion_vector["00000110"]=7;
    motion_vector["0000010110"]=8;
    motion_vector["0000010100"]=9;
    motion_vector["0000010010"]=10;
    motion_vector["00000100010"]=11;
    motion_vector["00000100000"]=12;
    motion_vector["00000011110"]=13;
    motion_vector["00000011100"]=14;
    motion_vector["00000011010"]=15;
    motion_vector["00000011000"]=16;
}

void build_dct_luminance()
{
    dct_luminance["100"]=0;
    dct_luminance["00"]=1;
    dct_luminance["01"]=2;
    dct_luminance["101"]=3;
    dct_luminance["110"]=4;
    dct_luminance["1110"]=5;
    dct_luminance["11110"]=6;
    dct_luminance["111110"]=7;
    dct_luminance["1111110"]=8;
}

void build_dct_chrominance()
{
    dct_chrominance["00"]=0;
    dct_chrominance["01"]=1;
    dct_chrominance["10"]=2;
    dct_chrominance["110"]=3;
    dct_chrominance["1110"]=4;
    dct_chrominance["11110"]=5;
    dct_chrominance["111110"]=6;
    dct_chrominance["1111110"]=7;
    dct_chrominance["11111110"]=8;
}

void build_dct_coeff()
{
    //加上"01"為error
    //run
    dct_coeff_run["10"]=-1;//end
    dct_coeff_run["1"]=0;//0,but only for coeff first
    dct_coeff_run["11"]=0;
    dct_coeff_run["011"]=1;
    dct_coeff_run["0100"]=0;
    dct_coeff_run["0101"]=2;
    dct_coeff_run["00101"]=0;
    dct_coeff_run["00111"]=3;
    dct_coeff_run["00110"]=4;
    dct_coeff_run["000110"]=1;
    dct_coeff_run["000111"]=5;
    dct_coeff_run["000101"]=6;
    dct_coeff_run["000100"]=7;
    dct_coeff_run["0000110"]=0;
    dct_coeff_run["0000100"]=2;
    dct_coeff_run["0000111"]=8;
    dct_coeff_run["0000101"]=9;
    dct_coeff_run["000001"]=-2;//escape;
    dct_coeff_run["00100110"]=0;
    dct_coeff_run["00100001"]=0;
    dct_coeff_run["00100101"]=1;
    dct_coeff_run["00100100"]=3;
    dct_coeff_run["00100111"]=10;
    dct_coeff_run["00100011"]=11;
    dct_coeff_run["00100010"]=12;
    dct_coeff_run["00100000"]=13;
    dct_coeff_run["0000001010"]=0;
    dct_coeff_run["0000001100"]=1;
    dct_coeff_run["0000001011"]=2;
    dct_coeff_run["0000001111"]=4;
    dct_coeff_run["0000001001"]=5;
    dct_coeff_run["0000001110"]=14;
    dct_coeff_run["0000001101"]=15;
    dct_coeff_run["0000001000"]=16;
    //page1 done
    dct_coeff_run["000000011101"]=0;
    dct_coeff_run["000000011000"]=0;
    dct_coeff_run["000000010011"]=0;
    dct_coeff_run["000000010000"]=0;
    dct_coeff_run["000000011011"]=1;
    dct_coeff_run["000000010100"]=2;
    dct_coeff_run["000000011100"]=3;
    dct_coeff_run["000000010010"]=4;
    dct_coeff_run["000000011110"]=6;
    dct_coeff_run["000000010101"]=7;
    dct_coeff_run["000000010001"]=8;
    dct_coeff_run["000000011111"]=17;
    dct_coeff_run["000000011010"]=18;
    dct_coeff_run["000000011001"]=19;
    dct_coeff_run["000000010111"]=20;
    dct_coeff_run["000000010110"]=21;
    dct_coeff_run["0000000011010"]=0;
    dct_coeff_run["0000000011001"]=0;
    dct_coeff_run["0000000011000"]=0;
    dct_coeff_run["0000000010111"]=0;
    dct_coeff_run["0000000010110"]=1;
    dct_coeff_run["0000000010101"]=1;
    dct_coeff_run["0000000010100"]=2;
    dct_coeff_run["0000000010011"]=3;
    dct_coeff_run["0000000010010"]=5;
    dct_coeff_run["0000000010001"]=9;
    dct_coeff_run["0000000010000"]=10;
    dct_coeff_run["0000000011111"]=22;
    dct_coeff_run["0000000011110"]=23;
    dct_coeff_run["0000000011101"]=24;
    dct_coeff_run["0000000011100"]=25;
    dct_coeff_run["0000000011011"]=26;
    //page2 done
    dct_coeff_run["00000000011111"]=0;
    dct_coeff_run["00000000011110"]=0;
    dct_coeff_run["00000000011101"]=0;
    dct_coeff_run["00000000011100"]=0;
    dct_coeff_run["00000000011011"]=0;
    dct_coeff_run["00000000011010"]=0;
    dct_coeff_run["00000000011001"]=0;
    dct_coeff_run["00000000011000"]=0;
    dct_coeff_run["00000000010111"]=0;
    dct_coeff_run["00000000010110"]=0;
    dct_coeff_run["00000000010101"]=0;
    dct_coeff_run["00000000010100"]=0;
    dct_coeff_run["00000000010011"]=0;
    dct_coeff_run["00000000010010"]=0;
    dct_coeff_run["00000000010001"]=0;
    dct_coeff_run["00000000010000"]=0;
    dct_coeff_run["000000000011000"]=0;
    dct_coeff_run["000000000010111"]=0;
    dct_coeff_run["000000000010110"]=0;
    dct_coeff_run["000000000010101"]=0;
    dct_coeff_run["000000000010100"]=0;
    dct_coeff_run["000000000010011"]=0;
    dct_coeff_run["000000000010010"]=0;
    dct_coeff_run["000000000010001"]=0;
    dct_coeff_run["000000000010000"]=0;
    dct_coeff_run["000000000011111"]=1;
    dct_coeff_run["000000000011110"]=1;
    dct_coeff_run["000000000011101"]=1;
    dct_coeff_run["000000000011100"]=1;
    dct_coeff_run["000000000011011"]=1;
    dct_coeff_run["000000000011010"]=1;
    dct_coeff_run["000000000011001"]=1;
    //page3 done
    dct_coeff_run["0000000000010011"]=1;
    dct_coeff_run["0000000000010010"]=1;
    dct_coeff_run["0000000000010001"]=1;
    dct_coeff_run["0000000000010000"]=1;
    dct_coeff_run["0000000000010100"]=6;
    dct_coeff_run["0000000000011010"]=11;
    dct_coeff_run["0000000000011001"]=12;
    dct_coeff_run["0000000000011000"]=13;
    dct_coeff_run["0000000000010111"]=14;
    dct_coeff_run["0000000000010110"]=15;
    dct_coeff_run["0000000000010101"]=16;
    dct_coeff_run["0000000000011111"]=27;
    dct_coeff_run["0000000000011110"]=28;
    dct_coeff_run["0000000000011101"]=29;
    dct_coeff_run["0000000000011100"]=30;
    dct_coeff_run["0000000000011011"]=31;
    ///////
    //level
    ///////
    dct_coeff_level["10"]=-1;//end
    dct_coeff_level["1"]=1;//only for coe first
    dct_coeff_level["11"]=1;
    dct_coeff_level["011"]=1;
    dct_coeff_level["0100"]=2;
    dct_coeff_level["0101"]=1;
    dct_coeff_level["00101"]=3;
    dct_coeff_level["00111"]=1;
    dct_coeff_level["00110"]=1;
    dct_coeff_level["000110"]=2;
    dct_coeff_level["000111"]=1;
    dct_coeff_level["000101"]=1;
    dct_coeff_level["000100"]=1;
    dct_coeff_level["0000110"]=4;
    dct_coeff_level["0000100"]=2;
    dct_coeff_level["0000111"]=1;
    dct_coeff_level["0000101"]=1;
    dct_coeff_level["000001"]=-2;//escape
    dct_coeff_level["00100110"]=5;
    dct_coeff_level["00100001"]=6;
    dct_coeff_level["00100101"]=3;
    dct_coeff_level["00100100"]=2;
    dct_coeff_level["00100111"]=1;
    dct_coeff_level["00100011"]=1;
    dct_coeff_level["00100010"]=1;
    dct_coeff_level["00100000"]=1;
    dct_coeff_level["0000001010"]=7;
    dct_coeff_level["0000001100"]=4;
    dct_coeff_level["0000001011"]=3;
    dct_coeff_level["0000001111"]=2;
    dct_coeff_level["0000001001"]=2;
    dct_coeff_level["0000001110"]=1;
    dct_coeff_level["0000001101"]=1;
    dct_coeff_level["0000001000"]=1;
    //page1 done
    dct_coeff_level["000000011101"]=8;
    dct_coeff_level["000000011000"]=9;
    dct_coeff_level["000000010011"]=10;
    dct_coeff_level["000000010000"]=11;
    dct_coeff_level["000000011011"]=5;
    dct_coeff_level["000000010100"]=4;
    dct_coeff_level["000000011100"]=3;
    dct_coeff_level["000000010010"]=3;
    dct_coeff_level["000000011110"]=2;
    dct_coeff_level["000000010101"]=2;
    dct_coeff_level["000000010001"]=2;
    dct_coeff_level["000000011111"]=1;
    dct_coeff_level["000000011010"]=1;
    dct_coeff_level["000000011001"]=1;
    dct_coeff_level["000000010111"]=1;
    dct_coeff_level["000000010110"]=1;
    dct_coeff_level["0000000011010"]=12;
    dct_coeff_level["0000000011001"]=13;
    dct_coeff_level["0000000011000"]=14;
    dct_coeff_level["0000000010111"]=15;
    dct_coeff_level["0000000010110"]=6;
    dct_coeff_level["0000000010101"]=7;
    dct_coeff_level["0000000010100"]=5;
    dct_coeff_level["0000000010011"]=4;
    dct_coeff_level["0000000010010"]=3;
    dct_coeff_level["0000000010001"]=2;
    dct_coeff_level["0000000010000"]=2;
    dct_coeff_level["0000000011111"]=1;
    dct_coeff_level["0000000011110"]=1;
    dct_coeff_level["0000000011101"]=1;
    dct_coeff_level["0000000011100"]=1;
    dct_coeff_level["0000000011011"]=1;
    //page2 done
    dct_coeff_level["00000000011111"]=16;
    dct_coeff_level["00000000011110"]=17;
    dct_coeff_level["00000000011101"]=18;
    dct_coeff_level["00000000011100"]=19;
    dct_coeff_level["00000000011011"]=20;
    dct_coeff_level["00000000011010"]=21;
    dct_coeff_level["00000000011001"]=22;
    dct_coeff_level["00000000011000"]=23;
    dct_coeff_level["00000000010111"]=24;
    dct_coeff_level["00000000010110"]=25;
    dct_coeff_level["00000000010101"]=26;
    dct_coeff_level["00000000010100"]=27;
    dct_coeff_level["00000000010011"]=28;
    dct_coeff_level["00000000010010"]=29;
    dct_coeff_level["00000000010001"]=30;
    dct_coeff_level["00000000010000"]=31;
    dct_coeff_level["000000000011000"]=32;
    dct_coeff_level["000000000010111"]=33;
    dct_coeff_level["000000000010110"]=34;
    dct_coeff_level["000000000010101"]=35;
    dct_coeff_level["000000000010100"]=36;
    dct_coeff_level["000000000010011"]=37;
    dct_coeff_level["000000000010010"]=38;
    dct_coeff_level["000000000010001"]=39;
    dct_coeff_level["000000000010000"]=40;
    dct_coeff_level["000000000011111"]=8;
    dct_coeff_level["000000000011110"]=9;
    dct_coeff_level["000000000011101"]=10;
    dct_coeff_level["000000000011100"]=11;
    dct_coeff_level["000000000011011"]=12;
    dct_coeff_level["000000000011010"]=13;
    dct_coeff_level["000000000011001"]=14;
    //page3 done
    dct_coeff_level["0000000000010011"]=15;
    dct_coeff_level["0000000000010010"]=16;
    dct_coeff_level["0000000000010001"]=17;
    dct_coeff_level["0000000000010000"]=18;
    dct_coeff_level["0000000000010100"]=3;
    dct_coeff_level["0000000000011010"]=2;
    dct_coeff_level["0000000000011001"]=2;
    dct_coeff_level["0000000000011000"]=2;
    dct_coeff_level["0000000000010111"]=2;
    dct_coeff_level["0000000000010110"]=2;
    dct_coeff_level["0000000000010101"]=2;
    dct_coeff_level["0000000000011111"]=1;
    dct_coeff_level["0000000000011110"]=1;
    dct_coeff_level["0000000000011101"]=1;
    dct_coeff_level["0000000000011100"]=1;
    dct_coeff_level["0000000000011011"]=1;
    
}
int sign(int t){
    return ((t>0)-(t<0));
}

void huff_tree(long i,int ptr,map<string,int> dist,unsigned char *buffer,long &i2, int &ptr2,int &result)
{
    //用來搜尋定義的huffman tree 主要傳入當前的bits及所要搜尋的tree
    string temp="";
    map<string, int>::iterator iter;
    
    do {
        temp+= hex_chars[buffer[i]>>(7-ptr)&0x01];
        i+=(ptr+1)/8;
        ptr=(ptr+1)%8;
        iter = dist.find(temp.c_str());
        if(iter != dist.end()) {
            result=iter->second;
            break;
        }
    }while(1);
    i2=i;
    ptr2=ptr;
}

void dct_next_tree(long i,int ptr,map<string,int> dist,unsigned char *buffer,long &i2, int &ptr2,int &result)
{
    //用來搜尋定義的huffman tree 主要傳入當前的bits及所要搜尋的tree
    string temp="";
    map<string, int>::iterator iter;
    
    do {
        temp+= hex_chars[buffer[i]>>(7-ptr)&0x01];
        i+=(ptr+1)/8;
        ptr=(ptr+1)%8;
        if ( strcmp(temp.c_str(),"1")==0 )continue;
        iter = dist.find(temp.c_str());
        if(iter != dist.end()) {
            result=iter->second;
            break;
        }
    }while(1);
    i2=i;
    ptr2=ptr;
}

void MB_type(long i,int ptr,int ft,unsigned char *buffer,long &i2, int &ptr2,string &result)
{
    //用來找macroblock type
    map<string, string>::iterator itr;
    string temp="";
    do {
        temp += hex_chars[buffer[i]>>(7-ptr)&0x01];
        i+=(ptr+1)/8;
        ptr=(ptr+1)%8;
        
        if (ft==1) {     // find macroblock_type by frame type
            itr=macroblock_type_I.find(temp.c_str());
            if(itr != macroblock_type_I.end()) {
                result=itr->second;
                break;
            }
        }else if (ft==2){
            itr=macroblock_type_P.find(temp.c_str());
            if(itr != macroblock_type_P.end()) {
                result=itr->second;
                break;
            }
        }else if (ft==3){
            itr=macroblock_type_B.find(temp.c_str());
            if(itr != macroblock_type_B.end()) {
                result=itr->second;
                break;
            }
        }else if (ft==4){
            itr=macroblock_type_D.find(temp.c_str());
            if(itr != macroblock_type_D.end()) {
                result=itr->second;
                break;
            }
        }else
            printf("frame number error");
    }while(1);
    i2=i;
    ptr2=ptr;
}
int ref_pel(int j,int m,int n,int right_for,int down_for,int right_half_for,int down_half_for, int p_n,int addr)
{
    //用來找past/future所要參考的pixel值
    int y=(addr/20)*16+down_for+m;
    int x=(addr%20)*16+right_for+n;
    int ch = 0;
    if (j<4){
        ch=0;
        if (j==1)x=x+8;
        if (j==2)y=y+8;
        if (j==3){
            x+=8;
            y+=8;
        }
    }else if(j==4){
        ch=1;
        x+=n;
        y+=m;
    }else if(j==5){
        ch=2;
        x+=n;
        y+=m;
    }
    while(x>319);
    while(y>239);
    
    
    if (!right_half_for&&!down_half_for)
        return gop[p_n][y][x][ch];
    if (!right_half_for&&down_half_for)
        return (gop[p_n][y][x][ch]+gop[p_n][y+1][x][ch])/2;
    if (right_half_for&&!down_half_for)
        return (gop[p_n][y][x][ch]+gop[p_n][y][x+1][ch])/2;
    if (right_half_for&&down_half_for)
        return (gop[p_n][y][x][ch]+gop[p_n][y+1][x][ch]+gop[p_n][y][x+1][ch]+gop[p_n][y+1][x+1][ch])/4;
    
    return 0;
}


int mb_width;
int mb_height;
long sequence(unsigned char *buffer , long i)
{
    int x[6];
    char temp[6];
    i+=4;
    
    while(1) {
        
        x[0] = buffer[i]/16;
        x[1] = buffer[i]%16;
        x[2] = buffer[i+1]/16;
        x[3] = buffer[i+1]%16;
        x[4] = buffer[i+2]/16;
        x[5] = buffer[i+2]%16;
        for (int j=0;j<6;j++){
            temp[j] = hex_chars[x[j]];
        }
        if ( strncmp(temp, "000001", 6) ==0 ) {
            break;
        }else {//讀影片的長跟寬
            int hor = (buffer[i]<<4)+(buffer[i+1]/16);
            int ver = ((buffer[i+1]%16)<<8)+buffer[i+2];
            mb_width=(hor+15)/16;
            mb_height=(ver+15)/16;
            //int flag = buffer[i+6]/64>>2;
            i+=3;
            //printf("hor:%d,ver:%d , flag:%d\n",hor,ver,flag);
            //printf("buffer:%d",buffer[i+4]);
            break;
        }
    }
    
    return i;
}

long group_of_pic(unsigned char *buffer , long i)
{
    int x[6];
    char temp[6];
    i+=4;
    
    while(1) {
        
        x[0] = buffer[i]/16;
        x[1] = buffer[i]%16;
        x[2] = buffer[i+1]/16;
        x[3] = buffer[i+1]%16;
        x[4] = buffer[i+2]/16;
        x[5] = buffer[i+2]%16;
        for (int j=0;j<6;j++){
            temp[j] = hex_chars[x[j]];
        }
        if ( strncmp(temp, "000001", 6) ==0 ) {
            break;
        }else {
            //buffer[i]:start of Group pic bitstream
            break;
        }
    }
    
    return i;
}


long pic_header(unsigned char *buffer , long i)
{
    int x[6];
    char temp[6];
    i+=4;
    
    while(1) {
        
        x[0] = buffer[i]/16;
        x[1] = buffer[i]%16;
        x[2] = buffer[i+1]/16;
        x[3] = buffer[i+1]%16;
        x[4] = buffer[i+2]/16;
        x[5] = buffer[i+2]%16;
        for (int j=0;j<6;j++){
            temp[j] = hex_chars[x[j]];
        }
        if ( strncmp(temp, "000001", 6) ==0 ) {
            break;
        }else {
            tsn = (buffer[i]<<2)+buffer[i+1]/64; // temperal sequence number
            ft = (buffer[i+1]%64)>>3; // frame type
            i++;
            //讀出當前的frame type及當前frame在播放影片時的位置
            
            for (int b=0;b<240;b++){
                for(int c=0;c<320;c++){
                    for(int d=0;d<3;d++){
                        gop[tsn][b][c][d]=0;
                    }
                }
            }//刷新用來存picture的buffer
            
            
            
            //
            i+=2;
            int ptr=5;
            if (ft==2||ft==3) {//若為P or B frame，則須在讀出picture的forward vector
                int next_bit =buffer[i]>>(7-ptr) &0x01;
                full_pel_forward_vector=next_bit;
                while (full_pel_forward_vector==1) printf("full_pel_forward_vector is 1\n");
                i+=(ptr+1)/8;ptr=(ptr+1)%8;
                int forward_f_code=0;
                for (int k=0;k<3;k++){
                    next_bit=buffer[i]>>(7-ptr) &0x01;
                    i+=(ptr+1)/8;ptr=(ptr+1)%8;
                    forward_f_code=(forward_f_code<<1)+next_bit;
                }
                forward_r_size=forward_f_code-1;
                forward_f=1<<forward_r_size;
                //printf("forward_f:%d,forward_r_size:%d\n",forward_f,forward_r_size);
            }
            if (ft==3) {//若為B frame，則需讀出picture的backward vector
                int next_bit =buffer[i]>>(7-ptr) &0x01;
                full_pel_backward_vector=next_bit;
                while (full_pel_backward_vector==1) printf("full_pel_backward_vector is 1\n");
                i+=(ptr+1)/8;ptr=(ptr+1)%8;
                int backward_f_code=0;
                for (int k=0;k<3;k++){
                    next_bit=buffer[i]>>(7-ptr) &0x01;
                    i+=(ptr+1)/8;ptr=(ptr+1)%8;
                    backward_f_code=(backward_f_code<<1)+next_bit;
                }
                backward_r_size=backward_f_code-1;
                backward_f=1<<backward_r_size;
                //printf("backward_f:%d,backward_r_size:%d\n",backward_f,backward_r_size);
            }
            //printf("tsn:%d , ft:%d\n",tsn,ft);
            //buffer[i]:start of Group pic bitstream
            break;
        }
    }
    
    return i;
}

void Macroblock(unsigned char *buffer , long i,int ptr,long &new_i,int &new_ptr);

long slice(unsigned char *buffer , long i)
{
    int x[6];
    char temp[6];
    i+=4;
    long check=i;
    pre_mb_addr=((buffer[i-1]&255)-1)*mb_width-1;
    past_intra_addr=-2;
    dct_dc_y_past=1024;
    dct_dc_cb_past=1024;
    dct_dc_cr_past=1024;
    
    recon_right_for_prev=0;
    recon_down_for_prev=0;
    recon_right_back_prev=0;
    recon_down_back_prev=0;
    
    int ptr=0;
    //每次進入slice值，要初始化參數
    while(1) {
        
        x[0] = buffer[i]/16;
        x[1] = buffer[i]%16;
        x[2] = buffer[i+1]/16;
        x[3] = buffer[i+1]%16;
        x[4] = buffer[i+2]/16;
        x[5] = buffer[i+2]%16;
        for (int j=0;j<6;j++){
            temp[j] = hex_chars[x[j]];
        }
        if ( strncmp(temp, "000001", 6) ==0 ) {
            break;
        }else {
            
            //quantizer_scale
            quan = buffer[i]>>3 ;
            //printf("quantizer_scale:%d\n",quan);
            ptr=5;
            int bf_sz=0;
            //slice extra
            int next_bit=buffer[i]>>(7-ptr) &0x01;
            i+=(ptr+1)/8;
            ptr=(ptr+1)%8;
            int next_byte=0;
            int extra_slice[10];
            while (next_bit==1) {//若下個bits為1，則在讀之後的8個bits做為extra information
                next_byte=((buffer[i]<<ptr) + (buffer[i+1]>>(7-ptr)))&0xff;
                extra_slice[bf_sz]=next_byte;
                bf_sz++;
                i++;
                next_bit=buffer[i]>>(7-ptr) &0x01;
                ptr=(ptr+1)%8;
            }
            
            for (int k=0;k<bf_sz;k++){
                printf("%d:%d\n",k,extra_slice[k]);
            }
            
            long new_i;
            int new_ptr;
            //Macroblock
            
            if (ft!=3){//若不為B frame，則更新用來參考的past picture及future picture
                ref_old=ref_new;
                ref_new=tsn;
            }
            
            while(1) {
                
                
                Macroblock(buffer,i,ptr,new_i,new_ptr);
                i=new_i;
                ptr=new_ptr;
                
                int sum=0;
                for (int j=0;j<23;j++){//檢查是否後23個bits皆為0（for enter next layer）
                    int next_bit=buffer[new_i]>>(7-new_ptr) &0x01;
                    new_i+=(new_ptr+1)/8;
                    new_ptr=(new_ptr+1)%8;
                    sum=sum+next_bit;
                }
                if (sum==0){
                    //printf("buffer:%d,buffer+1:%d,buffer+2:%d,buffer+3:%d\n",buffer[i],buffer[i+1],buffer[i+2],buffer[i+3]);
                    //printf("buffer:%d,buffer+1:%d,buffer+2:%d,buffer+3:%d\n",buffer[i+4],buffer[i+5],buffer[i+6],buffer[i+7]);
                    
                    imgwrite();
                    //解完一張picture之後，輸出為output.bmp
                    
                    //printf("frame_num:%d,frame type:%d\n",frame_num,ft);
                    //printf("tsn:%d,ref_old:%d,ref_new:%d\n",tsn,ref_old,ref_new);
                    break;
                }
            }
            break;
        }
    }
    return check;
}

void Macroblock(unsigned char *buffer , long i,int ptr,long &new_i,int &new_ptr)
{
    //Macroblock
    ///
    long i2;
    int ptr2;
    int next_bit = 0;
    //macroblock_address_increment
    
    
    
    int addr_inc;
    huff_tree(i, ptr, macro_address,buffer, i2, ptr2, addr_inc);//到huffman tree中找addr_increment
    i=i2;ptr=ptr2;
    //stuffing & escape
    while (addr_inc==34) printf("stuffing");
    int escape=0;
    while (addr_inc==35) {//進行escape處理
        escape+=1;
        //printf("escape:%d\n",escape);
        huff_tree(i, ptr, macro_address,buffer, i2, ptr2, addr_inc);
        i=i2;ptr=ptr2;
    }
    //addr_inc=addr_inc;
    addr_inc=addr_inc+33*escape;
    //
    pre_mb_addr+=addr_inc;//計算出當前block的編號
    //macroblock_type
    string macroblock_type;
    MB_type(i, ptr, ft, buffer, i2, ptr2 , macroblock_type);//找Macroblock的type
    i=i2;ptr=ptr2;
    //printf("\n:+: macroblock_number:%d :+:\n",pre_mb_addr+1);
    //printf("addr_inc:%d\n",addr_inc);
    //printf("macroblock_type:%s\n",macroblock_type.c_str());
    
    //1 : macroblock_quant
    //2 : macroblock_motion_forward
    //3 : macroblock_motion_backward
    //4 : macroblock_pattern
    //5 : macroblock_intra
    
    if (macroblock_type.c_str()[0]=='1'){
        //macroblock_quant
        quan=0;
        for (int k=0;k<5;k++){
            next_bit=buffer[i]>>(7-ptr) &0x01;
            i+=(ptr+1)/8;ptr=(ptr+1)%8;
            quan=(quan<<1)+next_bit;
        }
        //printf("macroblock_quant:%d\n",quan);
    }//找新的quantan值
    //macroblock_motion_forward
    if (macroblock_type.c_str()[1]=='1'){
        
        huff_tree(i, ptr, motion_vector,buffer, i2, ptr2, motion_horizontal_forward_code);
        i=i2;ptr=ptr2;
        
        
        
        if ((forward_f!=1)&&(motion_horizontal_forward_code!=0)){
            motion_horizontal_forward_r=0;
            for (int k=0;k<forward_r_size;k++){
                next_bit=buffer[i]>>(7-ptr) &0x01;
                i+=(ptr+1)/8;ptr=(ptr+1)%8;
                motion_horizontal_forward_r=(motion_horizontal_forward_r<<1)+next_bit;
            }
        }
        
        //vertical
        huff_tree(i, ptr, motion_vector,buffer, i2, ptr2, motion_vertical_forward_code);
        i=i2;ptr=ptr2;
        
        if ((forward_f!=1)&&(motion_vertical_forward_code!=0)){
            motion_vertical_forward_r=0;
            for (int k=0;k<forward_r_size;k++){
                next_bit=buffer[i]>>(7-ptr) &0x01;
                i+=(ptr+1)/8;ptr=(ptr+1)%8;
                motion_vertical_forward_r=(motion_vertical_forward_r<<1)+next_bit;
            }
        }
    }
    //macroblock_motion_backward
    if (macroblock_type.c_str()[2]=='1'){
        
        huff_tree(i, ptr, motion_vector,buffer, i2, ptr2, motion_horizontal_backward_code);
        i=i2;ptr=ptr2;
        
        
        motion_horizontal_backward_r=0;
        if ((backward_f!=1)&&(motion_horizontal_backward_code!=0)){
            
            //printf("motion_horizontal_backward_code:%d\n",motion_horizontal_backward_code);
            for (int k=0;k<backward_r_size;k++){
                next_bit=buffer[i]>>(7-ptr) &0x01;
                i+=(ptr+1)/8;ptr=(ptr+1)%8;
                motion_horizontal_backward_r=(motion_horizontal_backward_r<<1)+next_bit;
            }
            //printf("motion_horizontal_backward_r:%d\n",motion_horizontal_backward_r);
        }
        
        //vertical
        huff_tree(i, ptr, motion_vector,buffer, i2, ptr2, motion_vertical_backward_code);
        i=i2;ptr=ptr2;
        motion_vertical_backward_r=0;
        if ((backward_f!=1)&&(motion_vertical_backward_code!=0)){
            
            //printf("motion_vertical_backward_code:%d\n",motion_vertical_backward_code);
            for (int k=0;k<backward_r_size;k++){
                next_bit=buffer[i]>>(7-ptr) &0x01;
                i+=(ptr+1)/8;ptr=(ptr+1)%8;
                motion_vertical_backward_r=(motion_vertical_backward_r<<1)+next_bit;
            }
            //printf("motion_vertical_backward_r:%d\n",motion_vertical_backward_r);
        }
    }
    int cpb=0;
    if (macroblock_type.c_str()[3]=='1'){
        //macroblock_pattern
        huff_tree(i, ptr, macroblock_pattern,buffer, i2, ptr2, cpb);
        i=i2;ptr=ptr2;
        //printf("cpb:%d\n",cpb);
        //while(1) printf("macroblock_pattern\n");
    }//計算cpb值
    
    
    //derive pattern_code
    int pattern_code[6];
    for (int j=0;j<6;j++){
        pattern_code[j]=0;
        if (cpb&(1<<(5-j)))pattern_code[j]=1;
        if (macroblock_type.c_str()[4]=='1')pattern_code[j]=1;
    }
    
    int dct_block[6][64]={0};
    
    //利用huffman tree算出block的值
    for (int j=0;j<6;j++) {
        ////
        //block
        ////
        //for dct_zz
        if(pattern_code[j]){
            int run=0;
            int level=0;
            int count=0;
            int dct_size=0;
            int dct_diff = 0;
            //printf("%d's block\n",j);
            if(macroblock_type.c_str()[4]=='1'){ //macroblock_intra==1
                if (j<4){
                    //dct_dc_size_luminance
                    huff_tree(i, ptr, dct_luminance, buffer, i2, ptr2, dct_size);
                }else{
                    //dct_dc_size_chrominance
                    huff_tree(i, ptr, dct_chrominance, buffer, i2, ptr2, dct_size);
                
                }
                
                
                
                i=i2;ptr=ptr2;
                for (int k=0;k<dct_size;k++){
                    next_bit=buffer[i]>>(7-ptr) &0x01;
                    i+=(ptr+1)/8;ptr=(ptr+1)%8;
                    dct_diff=(dct_diff<<1)+next_bit;
                }
                if (dct_size>0){
                    if (dct_diff&(1<<(dct_size-1))) dct_block[j][0]=dct_diff;
                    else dct_block[j][0]=(-1<<dct_size)|(dct_diff+1);
                }else{
                    dct_block[j][0]=0;
                }
            }
            else{//macroblock_intra==0
                
                //dct_coeff_first
                
                huff_tree(i, ptr, dct_coeff_run, buffer, i2, ptr2, run);
                huff_tree(i, ptr, dct_coeff_level, buffer, i2, ptr2, level);
                i=i2;ptr=ptr2;
                
                while (run==-1) printf("here is problem\n");
                if (run==-2) { //escape case
                    run=0;
                    level=0;
                    //run
                    for (int k=0;k<6;k++){
                        next_bit=buffer[i]>>(7-ptr) &0x01;
                        i+=(ptr+1)/8;
                        ptr=(ptr+1)%8;
                        run=(run<<1)+next_bit;
                    }
                    //level
                    for (int k=0;k<8;k++){
                        next_bit=buffer[i]>>(7-ptr) &0x01;
                        i+=(ptr+1)/8;
                        ptr=(ptr+1)%8;
                        level=(level<<1)+next_bit;
                    }
                    if (level==0){
                        for (int k=0;k<8;k++){
                            next_bit=buffer[i]>>(7-ptr) &0x01;
                            i+=(ptr+1)/8;
                            ptr=(ptr+1)%8;
                            level=(level<<1)+next_bit;
                        }
                    }else if(level==128){
                        level=0;
                        for (int k=0;k<8;k++){
                            next_bit=buffer[i]>>(7-ptr) &0x01;
                            i+=(ptr+1)/8;
                            ptr=(ptr+1)%8;
                            level=(level<<1)+next_bit;
                        }
                        level=(-256+level);
                    }else if (level>128) {
                        level=-((level^0xff)+1);
                    }
                }else{
                    next_bit=buffer[i]>>(7-ptr) &0x01;
                    i+=(ptr+1)/8;ptr=(ptr+1)%8;
                    if (next_bit==0) level=level;
                    if (next_bit==1) level=-level;
                }
                count=run;
                dct_block[j][count]=level;
                //printf("dct_coeff_first,run:%d,level:%d\n",run,level);
                
                
            }
            
            if (ft!=4){
                //dct_coeff_next
                //printf("j=%d\ndct_diff=%d\n",j,dct_block[j][0]);
                
                while(1){
                    
                    //read until get a block
                    dct_next_tree(i, ptr, dct_coeff_run, buffer, i2, ptr2, run);
                    dct_next_tree(i, ptr, dct_coeff_level, buffer, i2, ptr2, level);
                    i=i2;ptr=ptr2;
                    
                    
                    if (run==-1) break;
                    if (run==-2) { //escape case
                        run=0;
                        level=0;
                        //run
                        for (int k=0;k<6;k++){
                            next_bit=buffer[i]>>(7-ptr) &0x01;
                            i+=(ptr+1)/8;
                            ptr=(ptr+1)%8;
                            run=(run<<1)+next_bit;
                        }
                        //level
                        for (int k=0;k<8;k++){
                            next_bit=buffer[i]>>(7-ptr) &0x01;
                            i+=(ptr+1)/8;
                            ptr=(ptr+1)%8;
                            level=(level<<1)+next_bit;
                        }
                        if (level==0){
                            for (int k=0;k<8;k++){
                                next_bit=buffer[i]>>(7-ptr) &0x01;
                                i+=(ptr+1)/8;
                                ptr=(ptr+1)%8;
                                level=(level<<1)+next_bit;
                            }
                        }else if(level==128){
                            level=0;
                            for (int k=0;k<8;k++){
                                next_bit=buffer[i]>>(7-ptr) &0x01;
                                i+=(ptr+1)/8;
                                ptr=(ptr+1)%8;
                                level=(level<<1)+next_bit;
                            }
                            level=(-256+level);
                        }else if (level>128) {
                            level=-((level^0xff)+1);
                        }
                    }else{
                        next_bit=buffer[i]>>(7-ptr) &0x01;
                        i+=(ptr+1)/8;ptr=(ptr+1)%8;
                        if (next_bit==0) level=level;
                        if (next_bit==1) level=-level;
                    }
                    count=count+run+1;
                    dct_block[j][count]=level;
                    //printf("dct_coeff_next ,run:%d,level:%d\n",run,level);
                }
                
            }
        }
    }
    
    
    int dct_recon[6][8][8]={0};
    /////
    ////reconstruct idct blocks
    ////y
    ///
    
    //在特定條件時初始化參數
    if (macroblock_type.c_str()[4]=='0'){
        dct_dc_y_past=1024;
        dct_dc_cb_past=1024;
        dct_dc_cr_past=1024;
    }
    if (pre_mb_addr==0){
        recon_right_for_prev=0;
        recon_down_for_prev=0;
        recon_right_back_prev=0;
        recon_down_back_prev=0;
    }
    if (pre_mb_addr==0||addr_inc>1){
        dct_dc_y_past=1024;
        dct_dc_cb_past=1024;
        dct_dc_cr_past=1024;
    }
    
    //intra block及p/b block的block處理
    if (macroblock_type.c_str()[4]=='1'){ //I-frame
        
        recon_right_for_prev=0;
        recon_down_for_prev=0;
        recon_right_back_prev=0;
        recon_down_back_prev=0;
        
        for (int j=0;j<6;j++) {
            for (int m=0;m<8;m++){
                for (int n=0;n<8;n++){
                    int k=scan[m][n];
                    dct_recon[j][m][n]=(2*dct_block[j][k]*quan*intra_quan[m][n])/16;
                    if( (dct_recon[j][m][n]&1)==0 ) dct_recon[j][m][n]=dct_recon[j][m][n]-sign(dct_recon[j][m][n]);
                    if( dct_recon[j][m][n]>2047   ) dct_recon[j][m][n]=2047;
                    if( dct_recon[j][m][n]<-2048  ) dct_recon[j][m][n]=-2048;
                }
            }
            if (j==0) {
                //y
                dct_recon[0][0][0]=dct_block[0][0]*8;
                if ( (pre_mb_addr-past_intra_addr)>1 )
                    dct_recon[0][0][0]=128*8+dct_recon[0][0][0];
                else
                    dct_recon[0][0][0]=dct_dc_y_past+dct_recon[0][0][0];
                dct_dc_y_past=dct_recon[0][0][0];
            }else if (j<4){
                dct_recon[j][0][0]=dct_dc_y_past + dct_block[j][0]*8;
                dct_dc_y_past = dct_recon[j][0][0] ;
            }else if (j==4){
    //    //cb
                dct_recon[4][0][0]=dct_block[4][0]*8;
                if ( (pre_mb_addr-past_intra_addr)>1 )
                    dct_recon[4][0][0]=128*8+dct_recon[4][0][0];
                else
                    dct_recon[4][0][0]=dct_dc_cb_past+dct_recon[4][0][0];
                dct_dc_cb_past=dct_recon[4][0][0];
            }else if (j==5){
    //    //cr
                dct_recon[5][0][0]=dct_block[5][0]*8;
                if ( (pre_mb_addr-past_intra_addr)>1 )
                    dct_recon[5][0][0]=128*8+dct_recon[5][0][0];
                else
                    dct_recon[5][0][0]=dct_dc_cr_past+dct_recon[5][0][0];
                dct_dc_cr_past=dct_recon[5][0][0];
            }
            IDCT(dct_recon[j]);
            for (int m=0;m<8;m++){
                for (int n=0;n<8;n++){
                    if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                    if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                }
            }
        }
        
    }else if (ft==2){// P frame
        
        if ( addr_inc>1 || macroblock_type.c_str()[1]=='0'){
            recon_right_for_prev=0;
            recon_down_for_prev=0;
        }
        
        int complement_horizontal_forward_r;
        int complement_vertical_forward_r;
        int right_little;
        int right_big;
        int down_little;
        int down_big;
        
        int max=(16*forward_f)-1;
        int min=(-16*forward_f);
        int new_vector;
        int recon_right_for=0;
        int recon_down_for=0;
        
        int right_for;
        int down_for;
        int right_half_for;
        int down_half_for;
        
        //printf("p_right:%d,p_down:%d\n",recon_right_for_prev,recon_down_for_prev);
        
        
        
        if(macroblock_type.c_str()[1]=='1') {
            if (forward_f==1 || motion_horizontal_forward_code==0){
                complement_horizontal_forward_r=0;
            }else{
                complement_horizontal_forward_r=forward_f-1-motion_horizontal_forward_r;
            }
            if (forward_f==1 || motion_vertical_forward_code==0){
                complement_vertical_forward_r=0;
            }else{
                complement_vertical_forward_r=forward_f-1-motion_vertical_forward_r;
            }
            right_little=motion_horizontal_forward_code*forward_f;
            if(right_little==0){
                right_big=0;
            }else{
                if(right_little>0){
                    right_little=right_little-complement_horizontal_forward_r;
                    right_big=right_little-32*forward_f;
                }else{
                    right_little=right_little+complement_horizontal_forward_r;
                    right_big=right_little+32*forward_f;
                }
            }
            down_little=motion_vertical_forward_code*forward_f;
            if(down_little==0){
                down_big=0;
            }else{
                if(down_little>0){
                    down_little=down_little-complement_vertical_forward_r;
                    down_big=down_little-32*forward_f;
                }else{
                    down_little=down_little+complement_vertical_forward_r;
                    down_big=down_little+32*forward_f;
                }
            }
            
            new_vector=recon_right_for_prev+right_little;
            if(new_vector<=max&&new_vector>=min){
                recon_right_for=recon_right_for_prev+right_little;
            }else{
                recon_right_for=recon_right_for_prev+right_big;
            }
            recon_right_for_prev=recon_right_for;
            while(full_pel_forward_vector)recon_right_for=recon_right_for<<1; //if>>while
            
            new_vector=recon_down_for_prev+down_little;
            if(new_vector<=max&&new_vector>=min){
                recon_down_for=recon_down_for_prev+down_little;
            }else{
                recon_down_for=recon_down_for_prev+down_big;
            }
            recon_down_for_prev=recon_down_for;
            while(full_pel_forward_vector)recon_down_for=recon_down_for<<1; //if>>while
        }
        
        //printf("MV(%d,%d)\n",recon_down_for,recon_right_for);
        //printf("mb_addr:%d\n",pre_mb_addr);
        for (int j=0;j<6;j++){
            if (pattern_code[j]){
                for(int m=0;m<8;m++){
                    for(int n=0;n<8;n++){
                        int k=scan[m][n];
                        dct_recon[j][m][n]=(((2*dct_block[j][k])+sign(dct_block[j][k]))*quan*16)/16; // *16,16:non_intra_quant[m][n]
                        if ((dct_recon[j][m][n]&1)==0) dct_recon[j][m][n]=dct_recon[j][m][n]-sign(dct_recon[j][m][n]);
                        if (dct_recon[j][m][n]>2047) dct_recon[j][m][n]=2047;
                        if (dct_recon[j][m][n]<-2048)dct_recon[j][m][n]=-2048;
                        if (dct_block[j][k]==0)dct_recon[j][m][n]=0;
                        
                    }
                }
                IDCT(dct_recon[j]);
//                for (int m=0;m<8;m++){
//                    for (int n=0;n<8;n++){
//                        if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
//                        if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
//                    }
//                }
            }
            if (j<4){
                right_for=recon_right_for>>1;
                down_for=recon_down_for>>1;
                right_half_for=recon_right_for-2*right_for;
                down_half_for=recon_down_for-2*down_for;
            }else{
                right_for=(recon_right_for/2)>>1;
                down_for=(recon_down_for/2)>>1;
                right_half_for=recon_right_for/2 - 2*right_for;
                down_half_for=recon_down_for/2 - 2*down_for;
            }
            
            if (macroblock_type.c_str()[1]=='0'){
                recon_right_for=0;
                recon_down_for=0;
            }
            for (int m=0;m<8;m++){
                for (int n=0;n<8;n++){
                    int ref=ref_pel(j,m,n,right_for,down_for,right_half_for,down_half_for,ref_old,pre_mb_addr);
                    dct_recon[j][m][n]+=ref;
                    if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                    if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                }
            }
        }
    }else if (ft==3){// B frame
        
        
        //for
        int recon_right_for=0;
        int recon_down_for=0;
        //back
        int recon_right_back=0;
        int recon_down_back=0;
        
        
        int complement_horizontal_forward_r;
        int complement_vertical_forward_r;
        int complement_horizontal_backward_r;
        int complement_vertical_backward_r;
        
        //skipblock
        if (addr_inc>1 && macroblock_type.c_str()[4]=='0'){
            int recon_right_for=past_for_r;
            int recon_down_for=past_for_d;
            //back
            int recon_right_back=past_back_r;
            int recon_down_back=past_back_d;
            
            for (int i=pre_mb_addr-addr_inc+1;i<pre_mb_addr;i++){
                //printf("skip block:%d\n",i);
                int dct_recon[6][8][8]={0};
                for (int j=0;j<6;j++){
                    
                    int right_for = 0, right_half_for = 0;
                    int down_for = 0, down_half_for = 0;
                    int right_back, right_half_back;
                    int down_back, down_half_back;
                    
                    if (macroblock_type.c_str()[1]=='1'){
                        if (j<4){
                            right_for=recon_right_for>>1;
                            down_for=recon_down_for>>1;
                            right_half_for=recon_right_for-2*right_for;
                            down_half_for=recon_down_for-2*down_for;
                        }else{
                            right_for=(recon_right_for/2)>>1;
                            down_for=(recon_down_for/2)>>1;
                            right_half_for=recon_right_for/2 - 2*right_for;
                            down_half_for=recon_down_for/2 - 2*down_for;
                        }
                        
                        if (macroblock_type.c_str()[2]=='1'){
                            
                            if (j<4){
                                right_back=recon_right_back>>1;
                                down_back=recon_down_back>>1;
                                right_half_back=recon_right_back-2*right_back;
                                down_half_back=recon_down_back-2*down_back;
                            }else{
                                right_back=(recon_right_back/2)>>1;
                                down_back=(recon_down_back/2)>>1;
                                right_half_back=recon_right_back/2 - 2*right_back;
                                down_half_back=recon_down_back/2 - 2*down_back;
                            }
                            for (int m=0;m<8;m++){
                                for (int n=0;n<8;n++){
                                    int f_ref=ref_pel(j,m,n,right_for,down_for,right_half_for,down_half_for,ref_old,i);
                                    int b_ref=ref_pel(j,m,n,right_back,down_back,right_half_back,down_half_back,ref_new,i);
                                    dct_recon[j][m][n]+=(f_ref+b_ref)/2;
                                    if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                                    if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                                }
                            }
                        }else{//macroblock_type.c_str()[2]=='0'
                            for (int m=0;m<8;m++){
                                for (int n=0;n<8;n++){
                                    int f_ref=ref_pel(j,m,n,right_for,down_for,right_half_for,down_half_for,ref_old,i);
                                    dct_recon[j][m][n]+=f_ref;
                                    if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                                    if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                                }
                            }
                        }
                    }else { //macroblock_type.c_str()[1]=='0'
                        if (j<4){
                            right_back=recon_right_back>>1;
                            down_back=recon_down_back>>1;
                            right_half_back=recon_right_back-2*right_back;
                            down_half_back=recon_down_back-2*down_back;
                        }else{
                            right_back=(recon_right_back/2)>>1;
                            down_back=(recon_down_back/2)>>1;
                            right_half_back=recon_right_back/2 - 2*right_back;
                            down_half_back=recon_down_back/2 - 2*down_back;
                        }
                        for (int m=0;m<8;m++){
                            for (int n=0;n<8;n++){
                                int b_ref=ref_pel(j,m,n,right_back,down_back,right_half_back,down_half_back,ref_new,i);
                                dct_recon[j][m][n]+=b_ref;
                                if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                                if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                            }
                        }
                    }
                }
                
                for (int j=0;j<6;j++){
                    for (int m=0;m<8;m++){
                        for(int n=0;n<8;n++){
                            all_recon[i][j][m][n]=dct_recon[j][m][n];
                        }
                    }
                }
            }
        }
        
        
        if (macroblock_type.c_str()[1]=='0' && macroblock_type.c_str()[2]=='0'){
            recon_right_for=recon_right_for_prev;
            recon_down_for=recon_down_for_prev;
            recon_right_back=recon_right_back_prev;
            recon_down_back=recon_down_back_prev;
        }
        
        if (macroblock_type.c_str()[1]=='1'){
            int new_vector;
            int right_little;
            int right_big;
            int down_little;
            int down_big;
            int max=(16*forward_f)-1;
            int min=(-16*forward_f);
            
            if (forward_f==1 || motion_horizontal_forward_code==0){
                complement_horizontal_forward_r=0;
            }else{
                complement_horizontal_forward_r=forward_f-1-motion_horizontal_forward_r;
            }
            if (forward_f==1 || motion_vertical_forward_code==0){
                complement_vertical_forward_r=0;
            }else{
                complement_vertical_forward_r=forward_f-1-motion_vertical_forward_r;
            }
            right_little=motion_horizontal_forward_code*forward_f;
            if(right_little==0){
                right_big=0;
            }else{
                if(right_little>0){
                    right_little=right_little-complement_horizontal_forward_r;
                    right_big=right_little-32*forward_f;
                }else{
                    right_little=right_little+complement_horizontal_forward_r;
                    right_big=right_little+32*forward_f;
                }
            }
            down_little=motion_vertical_forward_code*forward_f;
            if(down_little==0){
                down_big=0;
            }else{
                if(down_little>0){
                    down_little=down_little-complement_vertical_forward_r;
                    down_big=down_little-32*forward_f;
                }else{
                    down_little=down_little+complement_vertical_forward_r;
                    down_big=down_little+32*forward_f;
                }
            }
            
            new_vector=recon_right_for_prev+right_little;
            if(new_vector<=max&&new_vector>=min){
                recon_right_for=recon_right_for_prev+right_little;
            }else{
                recon_right_for=recon_right_for_prev+right_big;
            }
            recon_right_for_prev=recon_right_for;
            while(full_pel_forward_vector)recon_right_for=recon_right_for<<1; //if>>while
            
            new_vector=recon_down_for_prev+down_little;
            if(new_vector<=max&&new_vector>=min){
                recon_down_for=recon_down_for_prev+down_little;
            }else{
                recon_down_for=recon_down_for_prev+down_big;
            }
            recon_down_for_prev=recon_down_for;
            while(full_pel_forward_vector)recon_down_for=recon_down_for<<1; //if>>while
        }
        
        if (macroblock_type.c_str()[2]=='1'){
            int new_vector;
            int right_little;
            int right_big;
            int down_little;
            int down_big;
            
            int max=(16*backward_f)-1;
            int min=(-16*backward_f);
            
            if (backward_f==1 || motion_horizontal_backward_code==0){
                complement_horizontal_backward_r=0;
            }else{
                complement_horizontal_backward_r=backward_f-1-motion_horizontal_backward_r;
            }
            if (backward_f==1 || motion_vertical_backward_code==0){
                complement_vertical_backward_r=0;
            }else{
                complement_vertical_backward_r=backward_f-1-motion_vertical_backward_r;
            }
            right_little=motion_horizontal_backward_code*backward_f;
            if(right_little==0){
                right_big=0;
            }else{
                if(right_little>0){
                    right_little=right_little-complement_horizontal_backward_r;
                    right_big=right_little-32*backward_f;
                }else{
                    right_little=right_little+complement_horizontal_backward_r;
                    right_big=right_little+32*backward_f;
                }
            }
            down_little=motion_vertical_backward_code*backward_f;
            if(down_little==0){
                down_big=0;
            }else{
                if(down_little>0){
                    down_little=down_little-complement_vertical_backward_r;
                    down_big=down_little-32*backward_f;
                }else{
                    down_little=down_little+complement_vertical_backward_r;
                    down_big=down_little+32*backward_f;
                }
            }
            
            new_vector=recon_right_back_prev+right_little;
            if(new_vector<=max&&new_vector>=min){
                recon_right_back=recon_right_back_prev+right_little;
            }else{
                recon_right_back=recon_right_back_prev+right_big;
            }
            recon_right_back_prev=recon_right_back;
            while(full_pel_backward_vector)recon_right_back=recon_right_back<<1; //if>>while
            
            new_vector=recon_down_back_prev+down_little;
            if(new_vector<=max&&new_vector>=min){
                recon_down_back=recon_down_back_prev+down_little;
            }else{
                recon_down_back=recon_down_back_prev+down_big;
            }
            recon_down_back_prev=recon_down_back;
            while(full_pel_backward_vector)recon_down_back=recon_down_back<<1; //if>>while
        }
        
        past_for_r=recon_right_for;
        past_for_d=recon_down_for;
        past_back_r=recon_right_back;
        past_back_d=recon_down_back;
        for (int j=0;j<6;j++){
            
            int right_for = 0, right_half_for = 0;
            int down_for = 0, down_half_for = 0;
            int right_back, right_half_back;
            int down_back, down_half_back;
            
            if (pattern_code[j]){
                for(int m=0;m<8;m++){
                    for(int n=0;n<8;n++){
                        int k=scan[m][n];
                        dct_recon[j][m][n]=(((2*dct_block[j][k])+sign(dct_block[j][k]))*quan*16)/16; // *16,16:non_intra_quant[m][n]
                        if ((dct_recon[j][m][n]&1)==0) dct_recon[j][m][n]=dct_recon[j][m][n]-sign(dct_recon[j][m][n]);
                        if (dct_recon[j][m][n]>2047) dct_recon[j][m][n]=2047;
                        if (dct_recon[j][m][n]<-2048)dct_recon[j][m][n]=-2048;
                        if (dct_block[j][k]==0)dct_recon[j][m][n]=0;
                        
                    }
                }
                IDCT(dct_recon[j]);
            }
            
            if (macroblock_type.c_str()[1]=='1'){
                if (j<4){
                    right_for=recon_right_for>>1;
                    down_for=recon_down_for>>1;
                    right_half_for=recon_right_for-2*right_for;
                    down_half_for=recon_down_for-2*down_for;
                }else{
                    right_for=(recon_right_for/2)>>1;
                    down_for=(recon_down_for/2)>>1;
                    right_half_for=recon_right_for/2 - 2*right_for;
                    down_half_for=recon_down_for/2 - 2*down_for;
                }
            
                if (macroblock_type.c_str()[2]=='1'){
                
                    if (j<4){
                        right_back=recon_right_back>>1;
                        down_back=recon_down_back>>1;
                        right_half_back=recon_right_back-2*right_back;
                        down_half_back=recon_down_back-2*down_back;
                    }else{
                        right_back=(recon_right_back/2)>>1;
                        down_back=(recon_down_back/2)>>1;
                        right_half_back=recon_right_back/2 - 2*right_back;
                        down_half_back=recon_down_back/2 - 2*down_back;
                    }
                    for (int m=0;m<8;m++){
                        for (int n=0;n<8;n++){
                            int f_ref=ref_pel(j,m,n,right_for,down_for,right_half_for,down_half_for,ref_old,pre_mb_addr);
                            int b_ref=ref_pel(j,m,n,right_back,down_back,right_half_back,down_half_back,ref_new,pre_mb_addr);
                            dct_recon[j][m][n]+=(f_ref+b_ref)/2;
                            if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                            if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                        }
                    }
                }else{//macroblock_type.c_str()[2]=='0'
                    for (int m=0;m<8;m++){
                        for (int n=0;n<8;n++){
                            int f_ref=ref_pel(j,m,n,right_for,down_for,right_half_for,down_half_for,ref_old,pre_mb_addr);
                            dct_recon[j][m][n]+=f_ref;
                            if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                            if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                        }
                    }
                }
            }else { //macroblock_type.c_str()[1]=='0'
                if (j<4){
                    right_back=recon_right_back>>1;
                    down_back=recon_down_back>>1;
                    right_half_back=recon_right_back-2*right_back;
                    down_half_back=recon_down_back-2*down_back;
                }else{
                    right_back=(recon_right_back/2)>>1;
                    down_back=(recon_down_back/2)>>1;
                    right_half_back=recon_right_back/2 - 2*right_back;
                    down_half_back=recon_down_back/2 - 2*down_back;
                }
                for (int m=0;m<8;m++){
                    for (int n=0;n<8;n++){
                        int b_ref=ref_pel(j,m,n,right_back,down_back,right_half_back,down_half_back,ref_new,pre_mb_addr);
                        dct_recon[j][m][n]+=b_ref;
                        if (dct_recon[j][m][n]>255) dct_recon[j][m][n]=255;
                        if (dct_recon[j][m][n]<0) dct_recon[j][m][n]=0;
                    }
                }
            }
        }
    }

    
    //將當前block儲存
    for (int j=0;j<6;j++){
        for (int m=0;m<8;m++){
            for(int n=0;n<8;n++){
                all_recon[pre_mb_addr][j][m][n]=dct_recon[j][m][n];
            }
        }
    }
    past_block_type=macroblock_type;
    past_intra_addr=pre_mb_addr;
    
    new_i=i;
    new_ptr=ptr;
    
    if (ft==4){
        //end of macroblock
    }
    //printf("ptr:%d\n",ptr);
}

int main(int argc, const char * argv[]) {
    // insert code here...
    
//    if (argc < 2){
//        printf("Need File Name\n");
//        return -1;
//    }
    
    
    FILE *video;
    long sz;
    
    video = fopen("ipb_all.m1v","rb");
    if (video == NULL) {
        printf("CANNOT OPEN VIDEO\n");
        return 1;
    } else {
        printf("File read successfully\n");
        fseek(video, 0L, SEEK_END);
        sz = ftell(video);
        printf("size is %ld\n",sz);
        
    }
    
    unsigned char buffer[sz];
    
    fseek(video, 0, SEEK_SET);
    fread(buffer,sizeof(buffer), 1, video);
    
    //dictinary build
    build_macro_address();
    build_macroblock_type();
    build_macroblock_pattern();
    build_motion_vector();
    build_dct_luminance();
    build_dct_chrominance();
    build_dct_coeff();
    
//    map<string, int>::iterator iter;
//    iter = macro_address.find("010");
//    if(iter == macro_address.end()) {
//        printf("NOT FOUND\n");
//    }else {
//        printf("%d\n",iter->second);
//    }
    
    long i = 0;
    int x[8];
    char temp[8];
    int part = 0;
    
    int first=1;
//    int parts,parte;
    while( i <sizeof(buffer)) {
        x[0] = buffer[i]/16;
        x[1] = buffer[i]%16;
        x[2] = buffer[i+1]/16;
        x[3] = buffer[i+1]%16;
        x[4] = buffer[i+2]/16;
        x[5] = buffer[i+2]%16;
        x[6] = buffer[i+3]/16;
        x[7] = buffer[i+3]%16;
        for (int j=0;j<8;j++){
            temp[j] = hex_chars[x[j]];
        }
        //將讀進來的byte轉成hex format
        
        
        if (part ==0 & strcmp(temp,"000001B3") ==0 ) { //enter sequnce
            printf("Sequence Layer at %ld\n",i);
            i = sequence(buffer,i);
            part=1;
        }
        
        else if (part <= 1 & strcmp(temp,"000001B8") ==0 ) {//enter group of pic
            //printf("Group of picture at %ld\n",i);
            
            if (first!=1){
                printf("Group done\n");
                group_write();
            }
            first--;
            //每個group結束之後，依照播放的順序將img存成bmp(group_write()function)
            
            //刷新用來存圖片的buffer
            
            
            i = group_of_pic(buffer,i);
            part=2;
        }
        
        else if ( strcmp(temp, "000001B7") ==0 ) { // end of video
            group_write();
            printf(":+: End of Video at %ld :+:\n",i);
        }

        else if (part <= 4 & strcmp(temp, "00000100") ==0 ) {//enter picture layer
            //printf("Picture header at %ld\n",i);
            frame_num++;
            i = pic_header(buffer,i);
            part=0;
        }

        else if (part <= 3 & strcmp(temp, "00000101") ==0 ) {//enter slice layer
            
            //printf("Slice at %ld\n",i);
            
            i = slice(buffer,i);
            part=0;
            //printf("\nslice done\n");
        }
        
        
        i++;
    }
    
    
    fclose(video);
    return 0;
}
