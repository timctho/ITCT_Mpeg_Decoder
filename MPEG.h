#include <string.h>
#include <math.h>


#define picture_start_code      0x00000100
#define slice_start_code        0x00000100
#define user_data_start_code    0x000001b2
#define sequence_header_code    0x000001b3
#define sequence_error_code     0x000001b4
#define extension_start_code    0x000001b5
#define sequence_end_code       0x000001b7
#define group_start_code        0x000001b8
#define PI 3.14159265
#define C2 1.84775907     
#define C4 1.41421356    
#define C6 0.765366865   

static HWND hConWnd;
void Display_Video(char*,HWND=0,int=0,int=0,int=0,int=0,int=0,int=0,int=0,int=0);
HWND WindowHandle(void);
int buffer_count=-1;
int start_frame, end_frame;
float speed_factor;
bool is_loop;


void error(char*);
void gotoNextCode();


void decodeSequence();
void seqHeader();
void GOP();
void picture();
void decodeSlice();
void macroblock();
void block(int);

int Sign(int);
void fist_lum();
void next_lum();
void chrom_cb();
void chrom_cr();
void IDCT2();
void setHbitmap();

int table_2B1(int);
int table_2B2a(int);
int table_2B2b(int);
int table_2B2c(int);
int table_2B3(int);
int table_2B4(int);
int escape(int);
int table_2B5(int,int);
int table_2B5a(int);
int table_2B5b(int);
int table_2B5c(int,int);
int table_2B5d(int);
int table_2B5e(int);
int table_2B5f(int);

/*
global var usage
*/
FILE *f;
int offset = 0;
long long int mpeg_size,mpeg_count=0;
unsigned char* mpgBuffer;

/*
sequence layer var
*/
int horizontial_size;
int vertical_size;
int pel_aspect_ratio;
int picture_rate;
int bit_rate;
int marker_bit;
int vbv_buffer_size;
int constrained_parameter_flag;
int load_intra_quantizer_matrix;
int load_non_intra_quantizer_matrix;
//int extension_start_code;
int sequence_extension_data;
int user_data;


int mb_width,mb_height;
double height_width;
double pictures_per_second;


int time_code;
int closed_gop;
int broken_link;



int temporal_reference;
int picture_coding_type;
int vbv_delay;

int full_pixel_forward_vector;
int forward_f_code;
int forward_r_size;
int forward_f=0;
int full_pixel_backward_vector;
int backward_f_code;
int backward_r_size;
int backward_f=0;
int extra_bit_picture;
int extra_information_picture;
int picture_extension_data;


/*
decodeSlice layer var
*/
int quantizer_scale;
int extra_bit_slice;
int extra_information_slice;
int slice_vertical_position;
int previous_macroblock_address;

/*
macroblock layer var
*/
int macroblock_stuffing;
int macroblock_escape;
int macroblock_address_increment;
int macroblock_address;
int mb_row,mb_column;

int macroblock_quant;
int macroblock_motion_forward;
int macroblock_motion_backward;
int macroblock_pattern;
int macroblock_intra;
int motion_horzontal_forward_code;
int motion_horzontal_forward_r;
int motion_vertical_forward_code;
int motion_vertical_forward_r;
int motion_horzontal_backward_code;
int motion_horzontal_backward_r;
int motion_vertical_backward_code;
int motion_vertical_backward_r;
int coded_block_pattern;
int end_of_macroblock;
int pattern_code[6];
int past_intra_address;

int dct_dc_size_luminance;
int dct_dc_size_chrominance;
int dct_dc_differential;
int dct_coeff_first;
int dct_coeff_next;
int dct_dc_y_past;
int dct_dc_cb_past;
int dct_dc_cr_past;
int dct_recon[8][8];
int accumulate_zz;
int run;
int level;
int sign;


int once = 0;
double DataBuf[8][8];
double C(int x){
	if (x==0)return 1/(sqrt((double)2));
	else return 1;
}




char outname[256];    	
double macroblock_unit[3][16][16];
int i_buffer[3][240][320];
void setMacroBlock(int);
void setPicBuffer(int,int);


unsigned char i_buffer2[120][240][320][4];
