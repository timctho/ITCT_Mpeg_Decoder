#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <thread>
#include <omp.h>
#include <tchar.h>
#include <string.h>
#include <windows.h> // WinAPI
#include <time.h>
#include <assert.h>
#include "MPEG.h"
#include "static_table.h"
using namespace std;


void gotoNextCode(){     
		if (offset == 0)mpeg_count-=1;//fseek(f,-1,SEEK_CUR);     
		else offset = 0;
}

int Sign(int number){
		if (number==0)return 0;
		else if (number>0)return 1;
		else return -1;
}


void fist_lum(){//puts("first");
		int m,n,i;
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				i=scan[m][n];
				dct_recon[m][n] = (2 * dct_zz[i] * quantizer_scale * intra_quantizer_matrix[m][n])/16;
				if ((dct_recon[m][n]&1)==0)
					dct_recon[m][n] = dct_recon[m][n] - Sign(dct_recon[m][n]);
				if (dct_recon[m][n]>2047)dct_recon[m][n]=2047;
				if (dct_recon[m][n]<-2048)dct_recon[m][n]=-2048;
			}
		}
		dct_recon[0][0] = dct_zz[0]*8;
		if ((macroblock_address - past_intra_address)>1)dct_recon[0][0] = 128*8 + dct_recon[0][0];
		else dct_recon[0][0] = dct_dc_y_past + dct_recon[0][0];
		dct_dc_y_past = dct_recon[0][0];
}

void next_lum(){//puts("next");
		int m,n,i;
		//printf("%d %d\n",dct_zz[0],dct_recon[0][0]);
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				i=scan[m][n];
				dct_recon[m][n] = (2 * dct_zz[i] * quantizer_scale * intra_quantizer_matrix[m][n])/16;
				if ((dct_recon[m][n]&1)==0)
					dct_recon[m][n] = dct_recon[m][n] - Sign(dct_recon[m][n]);
				if (dct_recon[m][n]>2047)dct_recon[m][n]=2047;
				if (dct_recon[m][n]<-2048)dct_recon[m][n]=-2048;
			}
		}
		//printf("%d %d\n",dct_dc_y_past,dct_recon[0][0]);
		dct_recon[0][0] = dct_dc_y_past + dct_zz[0] * 8;
		dct_dc_y_past = dct_recon[0][0];
}

void chrom_cb(){
		int m,n,i;
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				i=scan[m][n];
				dct_recon[m][n] = (2 * dct_zz[i] * quantizer_scale * intra_quantizer_matrix[m][n])/16;
				if ((dct_recon[m][n]&1)==0)
					dct_recon[m][n] = dct_recon[m][n] - Sign(dct_recon[m][n]);
				if (dct_recon[m][n]>2047)dct_recon[m][n]=2047;
				if (dct_recon[m][n]<-2048)dct_recon[m][n]=-2048;
			}
		}
		dct_recon[0][0] = dct_zz[0]*8;
		if ((macroblock_address - past_intra_address)>1)dct_recon[0][0] = 128*8 + dct_recon[0][0];
		else dct_recon[0][0] = dct_dc_cb_past + dct_recon[0][0];
		dct_dc_cb_past = dct_recon[0][0];
}

void chrom_cr(){
		int m,n,i;
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				i=scan[m][n];
				dct_recon[m][n] = (2 * dct_zz[i] * quantizer_scale * intra_quantizer_matrix[m][n])/16;
				if ((dct_recon[m][n]&1)==0)
					dct_recon[m][n] = dct_recon[m][n] - Sign(dct_recon[m][n]);
				if (dct_recon[m][n]>2047)dct_recon[m][n]=2047;
				if (dct_recon[m][n]<-2048)dct_recon[m][n]=-2048;
			}
		}
		dct_recon[0][0] = dct_zz[0]*8;
		if ((macroblock_address - past_intra_address)>1)dct_recon[0][0] = 128*8 + dct_recon[0][0];
		else dct_recon[0][0] = dct_dc_cr_past + dct_recon[0][0];
		dct_dc_cr_past = dct_recon[0][0];
}


void IDCT(double F[],double f[]){
		int i;
		//puts("Before:");
		//for (i=0;i<8;i++)printf("%lf ",F[i]);printf("\n");
		for (i=0;i<8;i++)F[i]*=(4*cos(PI*i/16)/C(i));
		double a[8]; 	
		a[0] = F[0]/16.0;
		a[1] = F[4]/8.0;
		a[2] = (F[2] - F[6])/8.0;
		a[3] = (F[2] + F[6])/8.0;
		a[4] = (F[5] - F[3])/8.0;
		double temp1 = (F[1]+F[7])/8.0;
		double temp2 = (F[3]+F[5])/8.0;
		a[5] = temp1 - temp2;
		a[6] = (F[1]-F[7])/8.0;
		a[7] = temp1 + temp2;
		double b[8];
		b[0] = a[0];
		b[1] = a[1];
		b[2] = a[2]*C4;
		b[3] = a[3];
		b[4] = -(a[4]*C2+a[6]*C6);
		b[5] = a[5]*C4;
		b[6] = (-a[4]*C6+a[6]*C2);
		b[7] = a[7];
		double n[8];
		double temp3 = b[6] - b[7];
		n[0] = temp3 - b[5];
		n[1] = b[0]-b[1];
		n[2] = b[2]-b[3];
		n[3] = b[0]+b[1];
		n[4] = temp3;
		n[5] = b[4];
		n[6] = b[3];
		n[7] = b[7];
		double m[8];
		m[0] = n[7];
		m[1] = n[0];
		m[2] = n[4];
		m[3] = n[1]+n[2];
		m[4] = n[3]+n[6];
		m[5] = n[1]-n[2];
		m[6] = n[3]-n[6];
		m[7] = n[5]-n[0];
		//double f[8];
		f[0] = m[4]+m[0];
		f[1] = m[3]+m[2];
		f[2] = m[5]-m[1];
		f[3] = m[6]-m[7];
		f[4] = m[6]+m[7];
		f[5] = m[5]+m[1];
		f[6] = m[3]-m[2];
		f[7] = m[4]-m[0];
		
}

void IDCT2(){
		int m,n,mycount;
		double dct_recon_buf[8][8];
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				dct_recon_buf[m][n] = dct_recon[m][n];
			}//puts("");
		}//puts("");
		for (m=0;m<8;m++){
			IDCT(dct_recon_buf[m],DataBuf[m]);
		}
		double F[8],f[8];
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				F[n] = DataBuf[n][m];
			}
			IDCT(F,f);
			for (n=0;n<8;n++){
				DataBuf[n][m] = f[n]+128;
			}
		}
		//puts("DataBuf:");
		for (m=0;m<8;m++){
			for (n=0;n<8;n++){
				//DataUnit[mycount][i][j][m][n]=
				DataBuf[m][n]-=128;
				//printf("%.0lf ",DataBuf[m][n]);
			}//puts("");
		}//puts("");
}


void setMacroBlock(int i){
		int m,n,y=(i/2)*8,x=(i%2)*8;
		//printf ("x=%d y=%d\n",x,y);
		if (i<4){
			for (m=0;m<8;m++){
				for (n=0;n<8;n++){
					macroblock_unit[0][m+y][n+x] = DataBuf[m][n];
				}
			}
		} else {
			for (m=0;m<16;m++){
				for (n=0;n<16;n++){
					macroblock_unit[i-3][m][n] = DataBuf[m/2][n/2];
				}
			}
		}
		/*for (m=0;m<16;m++){
		  for (n=0;n<16;n++){
		  printf("%.0f ",macroblock_unit[2][m][n]);
		  }puts("");
		  }puts("");
		  system("pause");*/

}


void setPicBuffer(int row,int col){
		int i,j,k,m;
		for (i=0;i<16;i++){
			for (j=0;j<16;j++){
				i_buffer[0][row*16+i][col*16+j] = (int)(macroblock_unit[0][i][j]*1.167+1.596*macroblock_unit[2][i][j]-222.921);
				i_buffer[1][row*16+i][col*16+j] = (int)(macroblock_unit[0][i][j]*1.167-0.39176*macroblock_unit[1][i][j]-0.81328*macroblock_unit[2][i][j]+135.576);
				i_buffer[2][row*16+i][col*16+j] = (int)(macroblock_unit[0][i][j]*1.167+2.017*macroblock_unit[1][i][j]-276.836);
				
				// clamp color
				for (m=0;m<3;m++){
					if (i_buffer[m][row*16+i][col*16+j]<0)
						i_buffer[m][row*16+i][col*16+j]=0;
					else if (i_buffer[m][row*16+i][col*16+j]>255)
						i_buffer[m][row*16+i][col*16+j]=255;
				}
				}
				
			}
}





void error (char* mes){
	puts(mes);
	system("pause");
	exit(0);
}

int pic_count=-1;



void seqHeader(){
	//unsigned int* test =(unsigned int*)malloc(sizeof(unsigned int));
	unsigned char a,b,c,d;
	unsigned int e;
	int i,j;
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	//printf("%x\n",e);
	if (e!=sequence_header_code)error("seqHeader error");

	//fscanf(f,"%c%c%c",&a,&b,&c);
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	horizontial_size = (a<<4)|(b>>4);
	vertical_size = ((b&0x0f)<<8)|(c);
	mb_width = (horizontial_size+15)/16;
	mb_height = (vertical_size+15)/16;

	//fscanf(f,"%c",&a);
	a=mpgBuffer[mpeg_count++];
	pel_aspect_ratio = (a)>>4;
	picture_rate = (a&0x0f);

	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];

	bit_rate = (a<<10)|(b<<2)|((c&0xc0)>>6);
	marker_bit = (c&0x20)>>5;
	vbv_buffer_size = ((c&0x1f)<<3)|((d&0xf8)>>3);
	constrained_parameter_flag = (d&0x04)>>2;
	load_intra_quantizer_matrix = (d&0x02)>>1;

	offset = 6;
	if (load_intra_quantizer_matrix){
		for (i=0;i<64;i++){
			e = 0;
			for (j=0;j<8;j++){
				e |= ((d>>(7-offset))&0x1)<<(7-j);
				offset = (offset+1)%8;
				if (offset==0)d=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&d);
			}
			intra_quantizer_matrix[zzy[i]][zzx[i]]=e;
		}
	}
	load_non_intra_quantizer_matrix = (d&0x01);
	if (load_non_intra_quantizer_matrix){
		for (i=0;i<64;i++){
			e = 0;
			for (j=0;j<8;j++){
				e |= ((d>>(7-offset))&0x1)<<(7-j);
				offset = (offset+1)%8;
				if (offset==0)d=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&d);
			}
			non_intra_quantizer_matrix[zzy[i]][zzx[i]]=e;
		}
	}
	gotoNextCode();
	//offset = 0;
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	mpeg_count-=4;
	//fseek(f,-4,SEEK_CUR);

	if (e == extension_start_code){
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		if (e!=extension_start_code)error("extension_start_code error");
		//fscanf(f,"%d",&a);
		//extension_start_code = a;
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c",&b,&c,&d);
		e = (b<<16)|(c<<8)|(d);
		mpeg_count-=3;
		//fseek(f,-3,SEEK_CUR);
		while (e!=0x00000100){
			sequence_extension_data = mpgBuffer[mpeg_count++];//fgetc(f);
			//fscanf(f,"%c%c%c",&b,&c,&d);
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			e = (b<<16)|(c<<8)|(d);
			//fseek(f,-3,SEEK_CUR);
			mpeg_count-=3;
		}
		gotoNextCode();
		//offset = 0;
	}
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	mpeg_count-=4;
	//fseek(f,-4,SEEK_CUR);
	if (e == user_data_start_code){
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c",&b,&c,&d);
		e = (b<<16)|(c<<8)|(d);
		mpeg_count-=3;
		//fseek(f,-3,SEEK_CUR);
		while (e!=0x00000100){
			user_data = b=mpgBuffer[mpeg_count++];//fgetc(f);
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			//fscanf(f,"%c%c%c",&b,&c,&d);
			e = (b<<16)|(c<<8)|(d);
			//fseek(f,-3,SEEK_CUR);
			mpeg_count-=3;
		}
		gotoNextCode();
		//offset = 0;
	}
	/*printf ("hs=%d\nvs=%d\npar=%x\npr=%d\nbr=%d\nmb=%d\nvbv=%d\ncpf=%d\nliqm=%d\nlniqm=%d\n",
	  horizontial_size,vertical_size,pel_aspect_ratio,picture_rate,
	  bit_rate,marker_bit,vbv_buffer_size,constrained_parameter_flag,
	  load_intra_quantizer_matrix,load_non_intra_quantizer_matrix);
	  system("pause");*/
}

int picture_count=0;

void GOP(){
	unsigned char a,b,c,d;
	unsigned int e;
	int i,j,k;

	for (k=0;k<3;k++)
		for (i=0;i<16;i++)
			for (j=0;j<16;j++)
				i_buffer[k][i][j]=0;


	picture_count=0;
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	//printf("%x %x %x %x",a,b,c,d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	if (e!=group_start_code)error("GOP error");

	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	time_code = (a<<17)|(b<<9)|(c<<1)|(d>>7);
	closed_gop = (d&0x40)>>6;
	broken_link = (d&0x20)>>5;
	offset = 3;
	gotoNextCode();
	//offset = 0;


	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	mpeg_count-=4;
	//fseek(f,-4,SEEK_CUR);

	if (e == extension_start_code){
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		if (e!=extension_start_code)error("extension_start_code error");
		//fscanf(f,"%d",&a);
		//extension_start_code = a;
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c",&b,&c,&d);
		e = (b<<16)|(c<<8)|(d);
		mpeg_count-=3;
		//fseek(f,-3,SEEK_CUR);
		while (e!=0x000001){
			sequence_extension_data = mpgBuffer[mpeg_count++];//fgetc(f);
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			//fscanf(f,"%c%c%c",&b,&c,&d);
			e = (b<<16)|(c<<8)|(d);
			//fseek(f,-3,SEEK_CUR);
			mpeg_count-=3;
		}
		gotoNextCode();
		//offset = 0;
	}
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	//fseek(f,-4,SEEK_CUR);
	mpeg_count-=4;
	if (e == user_data_start_code){
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		if (e!=user_data_start_code)error("user_data_start_code error");

		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c",&b,&c,&d);
		e = (b<<16)|(c<<8)|(d);
		//fseek(f,-3,SEEK_CUR);
		mpeg_count-=3;
		while (e!=0x000001){
			user_data = mpgBuffer[mpeg_count++];//fgetc(f);
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			//fscanf(f,"%c%c%c",&b,&c,&d);
			e = (b<<16)|(c<<8)|(d);
			mpeg_count-=3;
			//fseek(f,-3,SEEK_CUR);
		}
		gotoNextCode();
		//offset = 0;
	}
	do {

		//printf("%d\n",++pic_count);
		++pic_count;
		++buffer_count;
		picture();

		

		for(int i=0;i<240;i++){
			for(int j=0;j<320;j++){
			    i_buffer2[pic_count][i][j][2] = (unsigned char)i_buffer[0][i][j];
		        i_buffer2[pic_count][i][j][1] = (unsigned char)i_buffer[1][i][j];
			    i_buffer2[pic_count][i][j][0] = (unsigned char)i_buffer[2][i][j];
				i_buffer2[pic_count][i][j][3] = (unsigned char)1;
			
			}
		}

		//char temp[50];
        //sprintf(temp,"%d.bmp",pic_count);
		//BCX_Bitmap(NULL,hConWnd,1,320,0,0,0);
        //Sleep(33);



		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		//fseek(f,-4,SEEK_CUR);
		mpeg_count-=4;
	}while(e==picture_start_code);     
	//printf ("tc=%d\ncg=%d\nbl=%d\n",time_code,closed_gop,broken_link);
	//system("pause");
}


void picture(){
	unsigned char a,b,c,d;
	unsigned int e;
	int i;
	once=0;
	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	//printf("%x %x %x %x",a,b,c,d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	//printf("e=%x\n",e);
	if (e!=picture_start_code)error("picture_start_code error");

	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	temporal_reference = (a<<2)|(b>>6);
	picture_coding_type = (b&0x38)>>3;
	vbv_delay = ((b&0x07)<<13)|(c<<5)|(d>>3);
	offset = 5;
	//puts("IPB frame");
	if (picture_coding_type==2||picture_coding_type==3){//puts("\tPB frame");
		e = 0;
		for (i=0;i<4;i++){
			e |= ((d >> (7-offset))&0x1)<<(3-i);
			offset = (offset+1)%8;
			if (offset==0)d=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&d);
		}
		full_pixel_forward_vector = (e>>3)&0x1;
		forward_f_code = e&0x7;
		forward_f = 1 << (forward_r_size = forward_f_code-1);

	}
	if (picture_coding_type==3){
		e = 0;
		for (i=0;i<4;i++){
			e |= ((d >> (7-offset))&0x1)<<(3-i);
			offset = (offset+1)%8;
			if (offset==0)d=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&d);
		}
		full_pixel_backward_vector = (e>>3)&0x1;
		backward_f_code = e&0x7;
		backward_f = 1 << (backward_r_size = backward_f_code-1);
	}
	//e=0;
	e = ((d >> (7-offset))&0x1);
	//offset = (offset+1)%8;
	//if (offset==0)fscanf(f,"%c",&d);
	while (e==0x01){
		extra_information_picture = 0;
		for (i=0;i<9;i++){
			if (i==0){
				extra_bit_picture = (d>>(7-offset))&0x1;
			}else {
				extra_information_picture |= ((d>>(7-offset))&0x1)<<(9-i);
			}
			offset = (offset+1)%8;
			if (offset==0)d=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&d);
		}
	}
	extra_bit_picture = (d>>(7-offset))&0x1;
	offset = (offset+1)%8;
	gotoNextCode();
	//offset = 0;

	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	mpeg_count-=4;
	//fseek(f,-4,SEEK_CUR);

	if (e==extension_start_code){
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		if (e!=extension_start_code)error("extension_start_code error");

		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c",&b,&c,&d);
		e = (b<<16)|(c<<8)|(d);
		mpeg_count-=3;
		//fseek(f,-3,SEEK_CUR);
		while (e!=0x000001){
			picture_extension_data = mpgBuffer[mpeg_count++];//fgetc(f);
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			//fscanf(f,"%c%c%c",&b,&c,&d);
			e = (b<<16)|(c<<8)|(d);
			mpeg_count-=3;
			//fseek(f,-3,SEEK_CUR);
		}
		gotoNextCode();
		//offset = 0;
	}

	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	mpeg_count-=4;
	//fseek(f,-4,SEEK_CUR);
	if (e==user_data_start_code){
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		if (e!=user_data_start_code)error("user_data_start_code error");

		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c",&b,&c,&d);
		e = (b<<16)|(c<<8)|(d);
		mpeg_count-=3;
		//fseek(f,-3,SEEK_CUR);
		while (e!=0x000001){
			user_data = mpgBuffer[mpeg_count++];//fgetc(f);
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			//fscanf(f,"%c%c%c",&b,&c,&d);
			e = (b<<16)|(c<<8)|(d);
			mpeg_count-=3;
			//fseek(f,-3,SEEK_CUR);
		}
		gotoNextCode();
		//offset = 0;
	}
	//printf ("tr=%d\npct=%d\nvbv=%d\n", temporal_reference,picture_coding_type,vbv_delay);
	//system("pause");
	macroblock_address = -1;
	do {
		decodeSlice();
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		mpeg_count-=4;
		//fseek(f,-4,SEEK_CUR);
	}while(e>=0x00000101 && e<=0x000001af);
}


void decodeSlice(){
	unsigned char a,b,c,d;
	unsigned int e;
	int i;
	past_intra_address = -2;
	dct_dc_y_past = 1024;
	dct_dc_cb_past = 1024;
	dct_dc_cr_past = 1024;

	a=mpgBuffer[mpeg_count++];
	b=mpgBuffer[mpeg_count++];
	c=mpgBuffer[mpeg_count++];
	d=mpgBuffer[mpeg_count++];
	//fscanf(f,"%c%c%c%c",&a,&b,&c,&d);
	e = (a<<24)|(b<<16)|(c<<8)|(d);
	if (!(e>=0x00000101 && e<=0x000001af))error("slice_start_code error");

	slice_vertical_position = e&0xff;
	previous_macroblock_address = (slice_vertical_position-1)*mb_width-1;

	//fscanf(f,"%c",&a);
	a=mpgBuffer[mpeg_count++];
	quantizer_scale = a>>3;
	offset = 5;
	while (((a<<offset)>>7)==0x01){
		extra_information_slice = 0;
		for (i=0;i<9;i++){
			if (i==0){
				extra_bit_slice = (d>>(7-offset))&0x1;
			}else {
				extra_information_slice |= ((d>>(7-offset))&0x1)<<(9-i);
			}
			offset = (offset+1)%8;
			if (offset==0)d=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
	}
	extra_bit_slice = (d>>(7-offset))&0x1;
	offset = (offset+1)%8;
	
	int mycount=0;
	do {
		macroblock();
		mpeg_count-=1;
		
		a=mpgBuffer[mpeg_count++];
		e = 0;
		for (i=0;i<23;i++){
			e |= ((a>>(7-offset))&0x1)<<(22-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		int back_shift = -1;
		int temp = 23;
		while (temp>offset){
			back_shift--;
			temp-=8;
		}
		offset = ((offset + 24 ) - (23)) % 8;
		mpeg_count+=back_shift;//fseek(f,back_shift,SEEK_CUR);
		a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);


	}while(((e)&0x7fffff)!=0x000000);
	gotoNextCode();
	//if (offset == 0)fseek(f,-1,SEEK_CUR);
	//offset = 0;
	//system("pause");
}


void macroblock(){//puts("1");
	unsigned char a,b,c,d;
	unsigned int e;
	int i,j,k;
	macroblock_address_increment = 0;

	
	for (k=0;k<3;k++)
		for (i=0;i<16;i++)
			for (j=0;j<16;j++)
				macroblock_unit[k][i][j]=0;

	mpeg_count-=1;
	a=mpgBuffer[mpeg_count++];
	
	e = 0;
	for (i=0;i<11;i++){
		e |= ((a>>(7-offset))&0x1)<<(10-i);
		
		offset = (offset+1)%8;
		if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
	}//puts("3");

	while ((e&0x7ff)==0xf){
		macroblock_stuffing = 0;
		for (i=0;i<11;i++){
			macroblock_stuffing |= ((a>>(7-offset))&0x1)<<(10-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		e = macroblock_stuffing;
	}
	


	while ((e&0x7ff)==0x8){
		macroblock_escape = 0;
		for (i=0;i<11;i++){
			macroblock_escape |= ((a>>(7-offset))&0x1)<<(10-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		e = macroblock_escape;
		macroblock_address_increment += 33;
	}
	/*
	   macroblock address increment
	   */
	int mai_leng = table_2B1(e);
	int code_leng = mai_leng&0xf;
	int back_shift = -1;
	int temp = 11-code_leng,mai_tmp=0;
	while (temp>offset){
		back_shift--;
		temp-=8;
	}
	offset = ((offset + 16 ) - (11 - code_leng)) % 8;
	
	mpeg_count+=back_shift;
	a=mpgBuffer[mpeg_count++];
	//printf("a=%x of=%d\n",a,offset);
	macroblock_address_increment += mai_leng>>8;
	previous_macroblock_address = macroblock_address = previous_macroblock_address + macroblock_address_increment;
	mb_row = macroblock_address / mb_width;
	mb_column = macroblock_address % mb_width;


	/*
	   macroblock type
	   */
	int VLC_code_leng;
	e = 0;
	//printf("offset=%d\n",offset);
	switch (picture_coding_type){
		case 1:{ /* I frame*/
			       for (i=0;i<2;i++){
				       e |= ((a>>(7-offset))&0x1)<<(1-i);
				       offset = (offset+1)%8;
				       if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			       }
			       VLC_code_leng = table_2B2a(e);
			       temp = 2;
			       break;
		       }case 2:{ /* P frame*/
			       for (i=0;i<6;i++){
				       e |= ((a>>(7-offset))&0x1)<<(5-i);
				       offset = (offset+1)%8;
				       if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			       }
			       VLC_code_leng = table_2B2b(e);
			       temp = 6;
			       break;
		       }case 3:{ /* B frame*/
			       for (i=0;i<6;i++){
				       e |= ((a>>(7-offset))&0x1)<<(5-i);
				       offset = (offset+1)%8;
				       if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			       }
			       VLC_code_leng = table_2B2c(e);
			       //printf("e=%x a=%x vlc=%x\n",e,a,VLC_code_leng);
			       temp = 6;
			       break;
		       }
	}
	int temp2 = temp;
	code_leng = VLC_code_leng&0xf;
	temp -= code_leng;
	back_shift = -1;
	while (temp>offset){
		back_shift--;
		temp-=8;
	}
	offset = ((offset + 16 ) - (temp2 - code_leng)) % 8;
	//fseek(f,back_shift,SEEK_CUR);
	//fscanf(f,"%c",&a);
	mpeg_count+=back_shift;
	a=mpgBuffer[mpeg_count++];
	macroblock_quant = (VLC_code_leng>>12)&0x1;
	macroblock_motion_forward = (VLC_code_leng>>11)&0x1;
	macroblock_motion_backward = (VLC_code_leng>>10)&0x1;
	macroblock_pattern = (VLC_code_leng>>9)&0x1;
	macroblock_intra = (VLC_code_leng>>8)&0x1;

	
	if (macroblock_quant){
		quantizer_scale = 0;
		for (i=0;i<5;i++){
			quantizer_scale |= ((a>>(7-offset))&0x1)<<(4-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
	}

	/*
	   macroblock motion forward
	   */
	if (macroblock_motion_forward){
		e = 0;
		for (i=0;i<11;i++){
			e |= ((a>>(7-offset))&0x1)<<(10-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		motion_horzontal_forward_code = table_2B4(e);
		code_leng = motion_horzontal_forward_code&0xff;
		motion_horzontal_forward_code >>= 8;
		motion_horzontal_forward_code -= 16;
		back_shift = -1;
		temp = 11 - code_leng;
		while (temp>offset){
			back_shift--;
			temp-=8;
		}
		offset = ((offset + 16 ) - (11 - code_leng)) % 8;
		//fseek(f,back_shift,SEEK_CUR);
		//fscanf(f,"%c",&a);
		mpeg_count+=back_shift;
		a=mpgBuffer[mpeg_count++];

		if ((forward_f!=1) && (motion_horzontal_forward_code!=0)){
			motion_horzontal_forward_r = 0;
			for (i=0;i<forward_r_size;i++){
				motion_horzontal_forward_r |= ((a>>(7-offset))&0x1)<<(forward_r_size-1-i);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			}
		}

		e = 0;
		for (i=0;i<11;i++){
			e |= ((a>>(7-offset))&0x1)<<(10-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		motion_vertical_forward_code = table_2B4(e);
		code_leng = motion_vertical_forward_code&0xff;
		motion_vertical_forward_code >>= 8;
		motion_vertical_forward_code -= 16;
		back_shift = -1;
		temp = 11 - code_leng;
		while (temp>offset){
			back_shift--;
			temp-=8;
		}
		offset = ((offset + 16 ) - (11 - code_leng)) % 8;
		//fseek(f,back_shift,SEEK_CUR);
		//fscanf(f,"%c",&a);
		mpeg_count+=back_shift;
		a=mpgBuffer[mpeg_count++];

		if ((forward_f!=1) && (motion_vertical_forward_code!=0)){
			motion_vertical_forward_r = 0;
			for (i=0;i<forward_r_size;i++){
				motion_vertical_forward_r |= ((a>>(7-offset))&0x1)<<(forward_r_size-1-i);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			}
		}



	}

	/*
	   macroblock motion backward
	   */
	if (macroblock_motion_backward){
		e = 0;
		for (i=0;i<11;i++){
			e |= ((a>>(7-offset))&0x1)<<(10-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		motion_horzontal_backward_code = table_2B4(e);
		code_leng = motion_horzontal_backward_code&0xff;
		motion_horzontal_backward_code >>= 8;
		motion_horzontal_backward_code -= 16;
		back_shift = -1;
		temp = 11 - code_leng;
		while (temp>offset){
			back_shift--;
			temp-=8;
		}
		offset = ((offset + 16 ) - (11 - code_leng)) % 8;
		//fseek(f,back_shift,SEEK_CUR);
		//fscanf(f,"%c",&a);
		mpeg_count+=back_shift;
		a=mpgBuffer[mpeg_count++];

		if ((backward_f!=1) && (motion_horzontal_backward_code!=0)){
			motion_horzontal_backward_r = 0;
			for (i=0;i<backward_r_size;i++){
				motion_horzontal_backward_r |= ((a>>(7-offset))&0x1)<<(backward_r_size-1-i);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			}
		}

		e = 0;
		for (i=0;i<11;i++){
			e |= ((a<<offset)>>7)<<(10-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		//printf("e=%x a=%x\n",e,a);
		motion_vertical_backward_code = table_2B4(e);
		code_leng = motion_vertical_backward_code&0xff;
		motion_vertical_backward_code >>= 8;
		motion_vertical_backward_code -=16;
		back_shift = -1;
		temp = 11 - code_leng;
		while (temp>offset){
			back_shift--;
			temp-=8;
		}
		offset = ((offset + 16 ) - (11 - code_leng)) % 8;
		//fseek(f,back_shift,SEEK_CUR);
		//fscanf(f,"%c",&a);
		mpeg_count+=back_shift;
		a=mpgBuffer[mpeg_count++];

		if ((backward_f!=1) && (motion_vertical_backward_code!=0)){
			motion_vertical_backward_r = 0;
			for (i=0;i<backward_r_size;i++){
				motion_vertical_backward_r |= ((a>>(7-offset))&0x1)<<(backward_r_size-1-i);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			}
		}

	}
	/*
	   macroblock pattern
	   */

	coded_block_pattern = 0;
	if (macroblock_pattern){
		e = 0;
		for (i=0;i<9;i++){
			e |= ((a>>(7-offset))&0x1)<<(8-i);
			offset = (offset+1)%8;
			if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		}
		coded_block_pattern = table_2B3(e);
		code_leng = coded_block_pattern&0xff;
		coded_block_pattern = coded_block_pattern>>8;
		back_shift = -1;
		temp = 9 - code_leng;
		while (temp>offset){
			back_shift--;
			temp-=8;
		}
		//printf ("e(use cpb)=%x a=%x offset=%d length=%d\n",e,a,offset,code_leng);
		offset = ((offset + 16 ) - (9 - code_leng)) % 8;
		//fseek(f,back_shift,SEEK_CUR);
		//fscanf(f,"%c",&a);
		mpeg_count+=back_shift;
		a=mpgBuffer[mpeg_count++];
		//printf ("e(use cpb)=%x a=%x offset=%d length=%d\n",e,a,offset,code_leng);
	}
	for (i=0;i<6;i++){
		pattern_code[i] = 0;
		if (coded_block_pattern&(1<<(5-i)) || macroblock_intra)pattern_code[i] = 1;
	}

	for (i=0;i<6;i++)
		block(i);
	past_intra_address = macroblock_address;
	setPicBuffer(mb_row,mb_column);
	int ii,jj;
	/*printf("macroblock (i=%d):\n",i);
	  for (ii=0;ii<16;ii++){
	  for (jj=0;jj<16;jj++){
	  printf ("%d ",i_buffer[0][ii+mb_row*16][jj+mb_column*16]);
	  }puts("");
	  }puts("");
	  system("pause");
	  */

	//system("pause");
	if (picture_coding_type==4){
		end_of_macroblock = ((a>>(7-offset))&0x1);
		offset = (offset+1)%8;
		if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
		if (end_of_macroblock!=1)error("end of macroblock error");
	}
	/*printf("ms=%d\nme=%d\nmai=%d\nqs=%d\nmhfc=%d\nmhfr=%d\nmvfc=%d\nmvfr=%d\ncbp=%d\n",
	  macroblock_stuffing,macroblock_escape,macroblock_address_increment,quantizer_scale,motion_horzontal_forward_code,
	  motion_horzontal_forward_r,motion_vertical_forward_code,motion_vertical_forward_r,coded_block_pattern);
	  system("pause");*/
}


void block(int i){
	unsigned char a,b,c,d;
	unsigned int e;
	int j;
	int code_leng,back_shift,temp;
	for (j=0;j<64;j++)dct_zz[j]=0;

	//fseek(f,-1,SEEK_CUR);
	//fscanf(f,"%c",&a);
	mpeg_count-=1;
	a=mpgBuffer[mpeg_count++];
	if (pattern_code[i]){

		if (macroblock_intra){/* dct dc size luminance */
			accumulate_zz=0;
			if (i<4){ 
				dct_zz[0] = 0;
				e = 0;
				for (j=0;j<7;j++){
					e |= ((a>>(7-offset))&0x1)<<(6-j);
					offset = (offset+1)%8;
					if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
				}
				dct_dc_size_luminance = table_2B5a(e);
				code_leng = dct_dc_size_luminance&0xf;
				dct_dc_size_luminance = dct_dc_size_luminance>>8;

				back_shift = -1;
				temp = 7 - code_leng;
				while (temp>offset){
					back_shift--;
					temp-=8;
				}
				offset = ((offset + 16 ) - (7 - code_leng)) % 8;
				//fseek(f,back_shift,SEEK_CUR);
				//fscanf(f,"%c",&a);
				mpeg_count+=back_shift;
				a=mpgBuffer[mpeg_count++];

				if (dct_dc_size_luminance>0){ /* dct dc differential */
					e = 0;
					for (j=0;j<dct_dc_size_luminance;j++){
						e |= ((a>>(7-offset))&0x1)<<(dct_dc_size_luminance-1-j);
						offset = (offset+1)%8;
						if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
					}
					dct_dc_differential = e;
					//if (dct_dc_differential>0){//puts("ddd>0");
					if (dct_dc_differential&(1<<(dct_dc_size_luminance-1)))dct_zz[0]=dct_dc_differential;
					else dct_zz[0]=(-1<<dct_dc_size_luminance)|(dct_dc_differential+1);
					//}
				}

			}else { 
				dct_zz[0] = 0;
				e = 0;
				for (j=0;j<8;j++){
					e |= ((a>>(7-offset))&0x1)<<(7-j);
					offset = (offset+1)%8;
					if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
				}
				dct_dc_size_chrominance = table_2B5b(e);
				code_leng = dct_dc_size_chrominance&0xf;
				dct_dc_size_chrominance = dct_dc_size_chrominance>>8;

				back_shift = -1;
				temp = 8 - code_leng;
				while (temp>offset){
					back_shift--;
					temp-=8;
				}
				offset = ((offset + 16 ) - (8 - code_leng)) % 8;
				//fseek(f,back_shift,SEEK_CUR);
				//fscanf(f,"%c",&a);
				mpeg_count+=back_shift;
				a=mpgBuffer[mpeg_count++];

				if (dct_dc_size_chrominance>0){ /* dct dc differential */
					e = 0;
					for (j=0;j<dct_dc_size_chrominance;j++){
						e |= ((a>>(7-offset))&0x1)<<(dct_dc_size_chrominance-1-j);
						offset = (offset+1)%8;
						if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
					}
					dct_dc_differential = e;
					//if (dct_dc_differential>0){
					if (dct_dc_differential&(1<<(dct_dc_size_chrominance-1)))dct_zz[0]=dct_dc_differential;
					else dct_zz[0]=(-1<<dct_dc_size_chrominance)|(dct_dc_differential+1);
					//}
				}
			}
		}else {
			e = 0;
			for (j=0;j<28;j++){
				e |= ((a>>(7-offset))&0x1)<<(27-j);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			}
			int cache = dct_coeff_first = table_2B5(e,0);// 0 for dct coeff first
			code_leng = dct_coeff_first&0xff;
			//dct_coeff_first = dct_coeff_first>>8;

			run = dct_coeff_first >> 24;
			level = (dct_coeff_first>>16)&0xff;
			sign = (dct_coeff_first>>8)&0xff;
			if (sign)level*=(-1);
			dct_zz [accumulate_zz=run] = level;

			back_shift = -1;
			temp = 28 - code_leng;
			while (temp>offset){
				back_shift--;
				temp-=8;
			}
			offset = ((offset + 32 ) - (28 - code_leng)) % 8;
			//fseek(f,back_shift,SEEK_CUR);
			//fscanf(f,"%c",&a);
			mpeg_count+=back_shift;
			a=mpgBuffer[mpeg_count++];


		}

		if (picture_coding_type!=4){
			e = 0;
			for (j=0;j<2;j++){
				e |= ((a>>(7-offset))&0x1)<<(1-j);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
			}
			back_shift = -1;
			temp = 2;
			while (temp>offset){
				back_shift--;
				temp-=8;
			}
			offset = ((offset + 8 ) - 2) % 8;
			//fseek(f,back_shift,SEEK_CUR);
			//fscanf(f,"%c",&a);
			mpeg_count+=back_shift;
			a=mpgBuffer[mpeg_count++];

			while (e!=0x2){
				e = 0;
				for (j=0;j<28;j++){
					e |= ((a>>(7-offset))&0x1)<<(27-j);
					offset = (offset+1)%8;
					if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
				}
				int cache = dct_coeff_next = table_2B5(e,1);/* 1 for dct coeff next*/
				code_leng = dct_coeff_next&0xff;
				//dct_coeff_next = dct_coeff_next>>8;

				run = dct_coeff_next >> 24;
				level = (dct_coeff_next>>16)&0xff;
				sign = (dct_coeff_next>>8)&0xff;
				if (sign)level*=(-1);
				accumulate_zz += (run+1);
				dct_zz[accumulate_zz]=level;

				back_shift = -1;
				temp = 28 - code_leng;
				while (temp>offset){
					back_shift--;
					temp-=8;
				}
				offset = ((offset + 32 ) - (28 - code_leng)) % 8;
				//fseek(f,back_shift,SEEK_CUR);
				//fscanf(f,"%c",&a);
				mpeg_count+=back_shift;
				a=mpgBuffer[mpeg_count++];

				e = 0;
				for (j=0;j<2;j++){
					e |= ((a>>(7-offset))&0x1)<<(1-j);
					offset = (offset+1)%8;
					if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
				}
				back_shift = -1;
				temp = 2;
				while (temp>offset){
					back_shift--;
					temp-=8;
				}
				offset = ((offset + 8 ) - (2)) % 8;
				//fseek(f,back_shift,SEEK_CUR);
				//fscanf(f,"%c",&a);
				mpeg_count+=back_shift;
				a=mpgBuffer[mpeg_count++];
			}
			e = 0;
			for (j=0;j<2;j++){
				e |= ((a>>(7-offset))&0x1)<<(1-j);
				offset = (offset+1)%8;
				if (offset==0)a=mpgBuffer[mpeg_count++];//fscanf(f,"%c",&a);
				
			}
			if (e!=0x2)error("end of block error");
			switch (i){
				case 0 : fist_lum();
					 break;
				case 1 :
					 next_lum();
					 break;
				case 2:
					 next_lum();
					 break;
				case 3:
					 next_lum();
					 break;
				case 4 : chrom_cb();
					 break;
				case 5 : chrom_cr();
					 break;
			}
			int ii,jj;

			IDCT2();

			setMacroBlock(i);

			}
		}
		//system("pause");
	}




	

void decodeSequence(){
	unsigned char a,b,c,d;
	unsigned int e;
	
	do{
		seqHeader();
		do{
			GOP();
			
			a=mpgBuffer[mpeg_count++];
			b=mpgBuffer[mpeg_count++];
			c=mpgBuffer[mpeg_count++];
			d=mpgBuffer[mpeg_count++];
			e = (a<<24)|(b<<16)|(c<<8)|(d);
			mpeg_count-=4;
			
		}while(e==group_start_code);
		
		a=mpgBuffer[mpeg_count++];
		b=mpgBuffer[mpeg_count++];
		c=mpgBuffer[mpeg_count++];
		d=mpgBuffer[mpeg_count++];
		e = (a<<24)|(b<<16)|(c<<8)|(d);
		//fseek(f,-4,SEEK_CUR);
		mpeg_count-=4;
	}while(e==sequence_header_code);
	if (e==sequence_end_code)
		printf("Decoding Complete.");
}
		
	

	



	int table_2B1(int code){
	int ccode = code>>3;
	int least = code&0x7;
	if ( ccode >= 0x80 && ccode <= 0xff) return (1<<8)|1;
	if ( ccode >= 0x60 && ccode <= 0x7f) return (2<<8)|3;
	if ( ccode >= 0x40 && ccode <= 0x5f) return (3<<8)|3;
	if ( ccode >= 0x30 && ccode <= 0x3f) return (4<<8)|4;
	if ( ccode >= 0x20 && ccode <= 0x2f) return (5<<8)|4;
	if ( ccode >= 0x18 && ccode <= 0x1f) return (6<<8)|5;
	if ( ccode >= 0x10 && ccode <= 0x17) return (7<<8)|5;
	if ( ccode >= 0x0e && ccode <= 0x0f) return (8<<8)|7;
	if ( ccode >= 0x0c && ccode <= 0x0d) return (9<<8)|7;
	if ( ccode == 0x0b) return (10<<8)|8;
	if ( ccode == 0x0a) return (11<<8)|8;
	if ( ccode == 0x09) return (12<<8)|8;
	if ( ccode == 0x08) return (13<<8)|8;
	if ( ccode == 0x07) return (14<<8)|8;
	if ( ccode == 0x06) return (15<<8)|8;
	if ( ccode == 0x05)
	{
		if ( least >= 0x6 && least <= 0x7) return (16<<8)|10;
		if ( least >= 0x4 && least <= 0x5) return (17<<8)|10;
		if ( least >= 0x2 && least <= 0x3) return (18<<8)|10;
		if ( least >= 0x0 && least <= 0x1) return (19<<8)|10;				  
	}
	if ( ccode == 0x04)
	{
		if ( least >= 0x6 && least <= 0x7) return (20<<8)|10;
		if ( least >= 0x4 && least <= 0x5) return (21<<8)|10;
		if ( least == 0x3) return (22<<8)|10;
		if ( least == 0x2) return (23<<8)|10;
		if ( least == 0x1) return (24<<8)|10;
		if ( least == 0x0) return (25<<8)|10;			  
	}
	if ( ccode == 0x03)
	{
		if ( least == 0x7) return (26<<8)|10;
		if ( least == 0x6) return (27<<8)|10;
		if ( least == 0x5) return (28<<8)|10;
		if ( least == 0x4) return (29<<8)|10;
		if ( least == 0x3) return (30<<8)|10;
		if ( least == 0x2) return (31<<8)|10;
		if ( least == 0x1) return (32<<8)|10;
		if ( least == 0x0) return (33<<8)|10;				  
	}
	if ( ccode == 0x01)
	{
		if ( least == 0x1) return (macroblock_stuffing<<8)|11;
		if ( least == 0x0) return (macroblock_escape<<8)|11;				  
	}
}

int table_2B2a(int code){
	switch (code){
		case 3 : return (0x01<<8)|1;
		case 2 : return (0x01<<8)|1;
		case 1 : return (0x11<<8)|2;
	}
}

int table_2B2b(int code){
	if ( code >= 0x20 && code <= 0x3f) return (0x0a<<8)|1;
	if ( code >= 0x10 && code <= 0x1f) return (0x02<<8)|2;
	if ( code >= 0x08 && code <= 0x0f) return (0x08<<8)|3;
	if ( code >= 0x06 && code <= 0x07) return (0x01<<8)|5;
	if ( code >= 0x04 && code <= 0x05) return (0x1a<<8)|5;
	if ( code >= 0x02 && code <= 0x03) return (0x12<<8)|5;
	if ( code == 0x01 ) return (0x11<<8)|6;
}

int table_2B2c(int code){
	if ( code >= 0x20 && code <= 0x2f) return (0x0c<<8)|2;
	if ( code >= 0x30 && code <= 0x3f) return (0x0e<<8)|2;
	if ( code >= 0x10 && code <= 0x17) return (0x04<<8)|3;
	if ( code >= 0x18 && code <= 0x1f) return (0x01<<8)|5;
	if ( code >= 0x08 && code <= 0x0b) return (0x08<<8)|4;
	if ( code >= 0x0c && code <= 0x0f) return (0x0a<<8)|4; 
	if ( code >= 0x06 && code <= 0x07) return (0x01<<8)|5;
	if ( code >= 0x04 && code <= 0x05) return (0x1e<<8)|5;
	if ( code == 0x03 ) return (0x1a<<8)|6;
	if ( code == 0x02 ) return (0x16<<8)|6;
	if ( code == 0x01 ) return (0x11<<8)|6;
}

int table_2B4(int code){
	int ccode = code << 1;
	if ( ccode >= 0x800 && ccode <= 0xfff) return ((0+16)<<8)|1;
	if ( ccode >= 0x400 && ccode <= 0x5ff) return ((1+16)<<8)|3;
	if ( ccode >= 0x600 && ccode <= 0x7ff) return ((-1+16)<<8)|3;
	if ( ccode >= 0x200 && ccode <= 0x2ff) return ((2+16)<<8)|4;
	if ( ccode >= 0x300 && ccode <= 0x3ff) return ((-2+16)<<8)|4;
	if ( ccode >= 0x100 && ccode <= 0x17f) return ((3+16)<<8)|5;  
	if ( ccode >= 0x180 && ccode <= 0x1ff) return ((-3+16)<<8)|5;
	if ( ccode >= 0x0c0 && ccode <= 0x0df) return ((4+16)<<8)|7;
	if ( ccode >= 0x0e0 && ccode <= 0x0ff) return ((-4+16)<<8)|7;
	if ( ccode >= 0x0a0 && ccode <= 0x0af) return ((5+16)<<8)|8;
	if ( ccode >= 0x0b0 && ccode <= 0x0bf) return ((-5+16)<<8)|8;
	if ( ccode >= 0x080 && ccode <= 0x08f) return ((6+16)<<8)|8;
	if ( ccode >= 0x090 && ccode <= 0x09f) return ((-6+16)<<8)|8;
	if ( ccode >= 0x060 && ccode <= 0x06f) return ((7+16)<<8)|8;
	if ( ccode >= 0x070 && ccode <= 0x07f) return ((-7+16)<<8)|8;
	if ( ccode >= 0x058 && ccode <= 0x05b) return ((8+16)<<8)|10;
	if ( ccode >= 0x05c && ccode <= 0x05f) return ((-8+16)<<8)|10;
	if ( ccode >= 0x050 && ccode <= 0x053) return ((9+16)<<8)|10;
	if ( ccode >= 0x054 && ccode <= 0x057) return ((-9+16)<<8)|10;
	if ( ccode >= 0x048 && ccode <= 0x04b) return ((10+16)<<8)|10;
	if ( ccode >= 0x04c && ccode <= 0x04f) return ((-10+16)<<8)|10;
	if ( ccode >= 0x044 && ccode <= 0x045) return ((11+16)<<8)|11;
	if ( ccode >= 0x046 && ccode <= 0x047) return ((-11+16)<<8)|11;
	if ( ccode >= 0x040 && ccode <= 0x041) return ((12+16)<<8)|11;
	if ( ccode >= 0x042 && ccode <= 0x043) return ((-12+16)<<8)|11;
	if ( ccode >= 0x03c && ccode <= 0x03d) return ((13+16)<<8)|11;
	if ( ccode >= 0x03e && ccode <= 0x03f) return ((-13+16)<<8)|11;
	if ( ccode >= 0x038 && ccode <= 0x039) return ((14+16)<<8)|11;
	if ( ccode >= 0x03a && ccode <= 0x03b) return ((-14+16)<<8)|11;
	if ( ccode >= 0x034 && ccode <= 0x035) return ((15+16)<<8)|11;
	if ( ccode >= 0x036 && ccode <= 0x037) return ((-15+16)<<8)|11;
	if ( ccode >= 0x030 && ccode <= 0x031) return ((16+16)<<8)|11;
	if ( ccode >= 0x032 && ccode <= 0x033) return ((-16+16)<<8)|11;    
}

int table_2B3(int code){
	int least = code & 1;
	int ccode = code >>1;
	if ( ccode >= 0xe0 && ccode <= 0xff) return (60<<8)|3;
	if ( ccode >= 0xd0 && ccode <= 0xdf) return (4<<8)|4;
	if ( ccode >= 0xc0 && ccode <= 0xcf) return (8<<8)|4;
	if ( ccode >= 0xb0 && ccode <= 0xbf) return (16<<8)|4;
	if ( ccode >= 0xa0 && ccode <= 0xaf) return (32<<8)|4;
	if ( ccode >= 0x98 && ccode <= 0x9f) return (12<<8)|5;
	if ( ccode >= 0x90 && ccode <= 0x97) return (48<<8)|5;
	if ( ccode >= 0x88 && ccode <= 0x8f) return (20<<8)|5;
	if ( ccode >= 0x80 && ccode <= 0x87) return (40<<8)|5;
	if ( ccode >= 0x78 && ccode <= 0x7f) return (28<<8)|5;
	if ( ccode >= 0x70 && ccode <= 0x77) return (44<<8)|5;
	if ( ccode >= 0x68 && ccode <= 0x6f) return (52<<8)|5;
	if ( ccode >= 0x60 && ccode <= 0x67) return (56<<8)|5;
	if ( ccode >= 0x58 && ccode <= 0x5f) return (1<<8)|5;
	if ( ccode >= 0x50 && ccode <= 0x57) return (61<<8)|5;
	if ( ccode >= 0x48 && ccode <= 0x4f) return (2<<8)|5;
	if ( ccode >= 0x40 && ccode <= 0x47) return (62<<8)|5;
	if ( ccode >= 0x3c && ccode <= 0x3f) return (24<<8)|6;
	if ( ccode >= 0x38 && ccode <= 0x3b) return (36<<8)|6;
	if ( ccode >= 0x34 && ccode <= 0x37) return (3<<8)|6;
	if ( ccode >= 0x30 && ccode <= 0x33) return (63<<8)|6;
	if ( ccode >= 0x2e && ccode <= 0x2f) return (5<<8)|7;
	if ( ccode >= 0x2c && ccode <= 0x2d) return (9<<8)|7;
	if ( ccode >= 0x2a && ccode <= 0x2b) return (17<<8)|7;
	if ( ccode >= 0x28 && ccode <= 0x29) return (33<<8)|7;
	if ( ccode >= 0x26 && ccode <= 0x27) return (6<<8)|7;
	if ( ccode >= 0x24 && ccode <= 0x25) return (10<<8)|7;
	if ( ccode >= 0x22 && ccode <= 0x23) return (18<<8)|7;
	if ( ccode >= 0x20 && ccode <= 0x21) return (34<<8)|7;
	if ( ccode == 0x1f )return (7<<8)|8;
	if ( ccode == 0x1e )return (11<<8)|8;
	if ( ccode == 0x1d )return (19<<8)|8;
	if ( ccode == 0x1c )return (35<<8)|8;
	if ( ccode == 0x1b )return (13<<8)|8;
	if ( ccode == 0x1a )return (49<<8)|8;
	if ( ccode == 0x19 )return (21<<8)|8;
	if ( ccode == 0x18 )return (41<<8)|8;
	if ( ccode == 0x17 )return (14<<8)|8;
	if ( ccode == 0x16 )return (50<<8)|8;
	if ( ccode == 0x15 )return (22<<8)|8;
	if ( ccode == 0x14 )return (42<<8)|8;
	if ( ccode == 0x13 )return (15<<8)|8;
	if ( ccode == 0x12 )return (51<<8)|8;
	if ( ccode == 0x11 )return (23<<8)|8;
	if ( ccode == 0x10 )return (43<<8)|8;
	if ( ccode == 0x0f )return (25<<8)|8;
	if ( ccode == 0x0e )return (37<<8)|8;
	if ( ccode == 0x0d )return (26<<8)|8;
	if ( ccode == 0x0c )return (38<<8)|8;
	if ( ccode == 0x0b )return (29<<8)|8;
	if ( ccode == 0x0a )return (45<<8)|8;
	if ( ccode == 0x09 )return (53<<8)|8;
	if ( ccode == 0x08 )return (57<<8)|8;
	if ( ccode == 0x07 )return (30<<8)|8;
	if ( ccode == 0x06 )return (46<<8)|8;
	if ( ccode == 0x05 )return (54<<8)|8;
	if ( ccode == 0x04 )return (58<<8)|8;
	if ( ccode == 0x03 )
	{
		switch (least){
			case 1 : return (31<<8)|9;
			case 0 : return (47<<8)|9;
		}
	}
	if ( ccode == 0x02 ){
		switch (least){
			case 1 : return (55<<8)|9;
			case 0 : return (59<<8)|9;
		}
	}
	if ( ccode == 0x01 ){
		switch (least){
			case 1 : return (27<<8)|9;
			case 0 : return (39<<8)|9;
		}
	}
}

int table_2B5a(int code){
	if ( code >= 0x40 && code <= 0x4f) return (0<<8)|3;
	if ( code >= 0x00 && code <= 0x1f) return (1<<8)|2;
	if ( code >= 0x20 && code <= 0x3f) return (2<<8)|2;
	if ( code >= 0x50 && code <= 0x5f) return (3<<8)|3;
	if ( code >= 0x60 && code <= 0x6f) return (4<<8)|3;
	if ( code >= 0x70 && code <= 0x77) return (5<<8)|4;
	if ( code >= 0x78 && code <= 0x7b) return (6<<8)|5;
	if ( code >= 0x7c && code <= 0x7d) return (7<<8)|6;
	if ( code == 0x7e) return (8<<8)|7;
}

int table_2B5b(int code){
	if ( code >= 0x00 && code <= 0x3f) return (0<<8)|2;
	if ( code >= 0x40 && code <= 0x7f) return (1<<8)|2;
	if ( code >= 0x80 && code <= 0xbf) return (2<<8)|2;
	if ( code >= 0xc0 && code <= 0xdf) return (3<<8)|3;
	if ( code >= 0xe0 && code <= 0xef) return (4<<8)|4;
	if ( code >= 0xf0 && code <= 0xf7) return (5<<8)|5;
	if ( code >= 0xf8 && code <= 0xfb) return (6<<8)|6;
	if ( code >= 0xfc && code <= 0xfd) return (7<<8)|7;
	if ( code == 0xfe) return (8<<8)|8;
}

int table_2B5(int code,int mode){ /* mode=0 for dct_coeff_first, mode=1 for dct_coeff_next*/
	if (((code&0x0ff00000)>>20)>1)return table_2B5c(code,mode);
	else if (((code&0x0fff0000)>>16)>7)return table_2B5d(code);
	else if (((code&0x0fff0000)>>16)>1)return table_2B5e(code);
	else return table_2B5f(code);
}

int table_2B5c(int code,int mode){/*return type: run:8bits, level:8bits, sign:8bits, length:8bits*/
	int ccode = code>>17;
	int least = ccode&0x7; 
	ccode >>= 3;
	if ((ccode>>7)==1){
		switch (mode){
			case 0 : return (0<<24)|(1<<16)|(((ccode>>6)&0x01)<<8)|2;
			case 1 : {
					 switch ((ccode>>6)&0x01){
						 case 0 : return 2;
						 case 1 : return (0<<24)|(1<<16)|(((ccode>>5)&0x01)<<8)|3;
					 }
				 }
		}
	} 
	else {
		if ( ccode >= 0x60 && ccode <= 0x7f) return (1<<24)|(1<<16)|(((ccode>>4)&0x01)<<8)|4;
		if ( ccode >= 0x40 && ccode <= 0x4f) return (0<<24)|(2<<16)|(((ccode>>3)&0x01)<<8)|5;
		if ( ccode >= 0x50 && ccode <= 0x5f) return (2<<24)|(1<<16)|(((ccode>>3)&0x01)<<8)|5;
		if ( ccode >= 0x28 && ccode <= 0x2f) return (0<<24)|(3<<16)|(((ccode>>2)&0x01)<<8)|6;
		if ( ccode >= 0x38 && ccode <= 0x3f) return (3<<24)|(1<<16)|(((ccode>>2)&0x01)<<8)|6;
		if ( ccode >= 0x30 && ccode <= 0x37) return (4<<24)|(1<<16)|(((ccode>>2)&0x01)<<8)|6;
		if ( ccode >= 0x18 && ccode <= 0x1b) return (1<<24)|(2<<16)|(((ccode>>1)&0x01)<<8)|7;
		if ( ccode >= 0x1c && ccode <= 0x1f) return (5<<24)|(1<<16)|(((ccode>>1)&0x01)<<8)|7;
		if ( ccode >= 0x14 && ccode <= 0x17) return (6<<24)|(1<<16)|(((ccode>>1)&0x01)<<8)|7;
		if ( ccode >= 0x10 && ccode <= 0x13) return (7<<24)|(1<<16)|(((ccode>>1)&0x01)<<8)|7;
		if ( ccode >= 0x0c && ccode <= 0x0d) return (0<<24)|(4<<16)|(((ccode)&0x01)<<8)|8;
		if ( ccode >= 0x08 && ccode <= 0x09) return (2<<24)|(2<<16)|(((ccode)&0x01)<<8)|8;
		if ( ccode >= 0x0e && ccode <= 0x0f) return (8<<24)|(1<<16)|(((ccode)&0x01)<<8)|8;
		if ( ccode >= 0x0a && ccode <= 0x0b) return (9<<24)|(1<<16)|(((ccode)&0x01)<<8)|8;
		if ( ccode >= 0x04 && ccode <= 0x07) return escape(code);
		if ( ccode == 0x26 ) return (0<<24)|(5<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x21 ) return (0<<24)|(6<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x25 ) return (1<<24)|(3<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x24 ) return (3<<24)|(2<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x27 ) return (10<<24)|(1<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x23 ) return (11<<24)|(1<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x22 ) return (12<<24)|(1<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x20 ) return (13<<24)|(1<<16)|(((least>>2)&0x01)<<8)|9;
		if ( ccode == 0x02 ) 
		{
			if ( least >= 0x4 && least <= 0x5) return (0<<24)|(7<<16)|(((least)&0x01)<<8)|11;
			if ( least >= 0x6 && least <= 0x7) return (2<<24)|(3<<16)|(((least)&0x01)<<8)|11;
			if ( least >= 0x2 && least <= 0x3) return (5<<24)|(2<<16)|(((least)&0x01)<<8)|11;
			if ( least >= 0x0 && least <= 0x1) return (16<<24)|(1<<16)|(((least)&0x01)<<8)|11;
		}
		if ( ccode == 0x03 ) 
		{
			if ( least >= 0x4 && least <= 0x5) return (14<<24)|(1<<16)|(((least)&0x01)<<8)|11;
			if ( least >= 0x6 && least <= 0x7) return (4<<24)|(2<<16)|(((least)&0x01)<<8)|11;
			if ( least >= 0x2 && least <= 0x3) return (15<<24)|(1<<16)|(((least)&0x01)<<8)|11;
			if ( least >= 0x0 && least <= 0x1) return (1<<24)|(4<<16)|(((least)&0x01)<<8)|11;
		}
	}   
}

int table_2B5d(int code){/*return type: run:8bits, level:8bits, sign:8bits, length:8bits*/
	int ccode = code >> 14;
	int least = ccode & 0x3;
	ccode >>= 2;

	switch (ccode){
		case 0x01d : return (0<<24)|(8<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x018 : return (0<<24)|(9<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x013 : return (0<<24)|(10<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x010 : return (0<<24)|(11<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x01b : return (1<<24)|(5<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x014 : return (2<<24)|(4<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x01c : return (3<<24)|(3<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x012 : return (4<<24)|(3<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x01e : return (6<<24)|(2<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x015 : return (7<<24)|(2<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x011 : return (8<<24)|(2<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x01f : return (17<<24)|(1<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x01a : return (18<<24)|(1<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x019 : return (19<<24)|(1<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x017 : return (20<<24)|(1<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x016 : return (21<<24)|(1<<16)|(((least>>1)&0x1)<<8)|13;
		case 0x00d : {
				     switch (least>>1){
					     case 0: return (0<<24)|(12<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (26<<24)|(1<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x00c : {
				     switch (least>>1){
					     case 0: return (0<<24)|(14<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (0<<24)|(13<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x00b : {
				     switch (least>>1){
					     case 0: return (1<<24)|(6<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (0<<24)|(15<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x00a : {
				     switch (least>>1){
					     case 0: return (2<<24)|(5<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (1<<24)|(7<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x009 : {
				     switch (least>>1){
					     case 0: return (5<<24)|(3<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (3<<24)|(4<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x008 : {
				     switch (least>>1){
					     case 0: return (10<<24)|(2<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (9<<24)|(2<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x00f : {
				     switch (least>>1){
					     case 0: return (23<<24)|(1<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (22<<24)|(1<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
		case 0x00e : {
				     switch (least>>1){
					     case 0: return (25<<24)|(1<<16)|(((least)&0x1)<<8)|14;
					     case 1: return (24<<24)|(1<<16)|(((least)&0x1)<<8)|14;
				     }
			     }
	}
}

int table_2B5e(int code){/*return type: run:8bits, level:8bits, sign:8bits, length:8bits*/
	int ccode = code >> 12;
	if ( ccode >= 0x007c && ccode <= 0x007f) return (0<<24)|(16<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0078 && ccode <= 0x007b) return (0<<24)|(17<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0074 && ccode <= 0x0077) return (0<<24)|(18<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0070 && ccode <= 0x0073) return (0<<24)|(19<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x006c && ccode <= 0x006f) return (0<<24)|(20<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0068 && ccode <= 0x006b) return (0<<24)|(21<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0064 && ccode <= 0x0067) return (0<<24)|(22<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0060 && ccode <= 0x0063) return (0<<24)|(23<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x005c && ccode <= 0x005f) return (0<<24)|(24<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0058 && ccode <= 0x005b) return (0<<24)|(25<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0054 && ccode <= 0x0057) return (0<<24)|(26<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0050 && ccode <= 0x0053) return (0<<24)|(27<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x004c && ccode <= 0x004f) return (0<<24)|(28<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0048 && ccode <= 0x004b) return (0<<24)|(29<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0044 && ccode <= 0x0047) return (0<<24)|(30<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0040 && ccode <= 0x0043) return (0<<24)|(31<<16)|(((ccode>>1)&0x1)<<8)|15;
	if ( ccode >= 0x0030 && ccode <= 0x0031) return (0<<24)|(32<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x002e && ccode <= 0x002f) return (0<<24)|(33<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x002c && ccode <= 0x002d) return (0<<24)|(34<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x002a && ccode <= 0x002b) return (0<<24)|(35<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0028 && ccode <= 0x0029) return (0<<24)|(36<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0026 && ccode <= 0x0027) return (0<<24)|(37<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0024 && ccode <= 0x0025) return (0<<24)|(38<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0022 && ccode <= 0x0023) return (0<<24)|(39<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0020 && ccode <= 0x0021) return (0<<24)|(40<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x003e && ccode <= 0x003f) return (1<<24)|(8<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x003c && ccode <= 0x003d) return (1<<24)|(9<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x003a && ccode <= 0x003b) return (1<<24)|(10<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0038 && ccode <= 0x0039) return (1<<24)|(11<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0036 && ccode <= 0x0037) return (1<<24)|(12<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0034 && ccode <= 0x0035) return (1<<24)|(13<<16)|(((ccode)&0x1)<<8)|16;
	if ( ccode >= 0x0032 && ccode <= 0x0033) return (1<<24)|(14<<16)|(((ccode)&0x1)<<8)|16;
}

int table_2B5f(int code){/*return type: run:8bits, level:8bits, sign:8bits, length:8bits*/
	int ccode = code >> 12;
	int least = (code>>11)&0x1;

	switch (ccode){
		case 0x0013 : return (1<<24)|(15<<16)|(least<<8)|17;
		case 0x0012 : return (1<<24)|(16<<16)|(least<<8)|17;
		case 0x0011 : return (1<<24)|(17<<16)|(least<<8)|17;
		case 0x0010 : return (1<<24)|(18<<16)|(least<<8)|17;
		case 0x0014 : return (6<<24)|(3<<16)|(least<<8)|17;
		case 0x001a : return (11<<24)|(2<<16)|(least<<8)|17;
		case 0x0019 : return (12<<24)|(2<<16)|(least<<8)|17;
		case 0x0018 : return (13<<24)|(2<<16)|(least<<8)|17;
		case 0x0017 : return (14<<24)|(2<<16)|(least<<8)|17;
		case 0x0016 : return (15<<24)|(2<<16)|(least<<8)|17;
		case 0x0015 : return (16<<24)|(2<<16)|(least<<8)|17;
		case 0x001f : return (27<<24)|(1<<16)|(least<<8)|17;
		case 0x001e : return (28<<24)|(1<<16)|(least<<8)|17;
		case 0x001d : return (29<<24)|(1<<16)|(least<<8)|17;
		case 0x001c : return (30<<24)|(1<<16)|(least<<8)|17;
		case 0x001b : return (31<<24)|(1<<16)|(least<<8)|17;
	}
}

int escape(int code){/*return type: run:8bits, level:8bits, sign:8bits, length:8bits*/
	int run = (code >> 16)&0x3f;
	int ccode = code&0xffff;
	int sign = ccode>>15;
	int level,length=20;
	switch (sign){
		case 0 : {
				 if (((ccode>>8&0xff))>0){
					 level = (ccode>>8)&0xff;
			         break;
				 }else{
					 level = ccode&0xff;
					 length+=8;
			         break;
				 }
				 break;
			 }
		case 1:
			 if (((ccode>>8)&0xff)>0x80){
				 level = (((~(ccode>>8))+1)&0xff);
				 //printf("\t\t\t(a)level=%d\n",level);
				 break;
			 }else{
				 level = (((~ccode)+1)&0xff);
				 //printf("\t\t\t(b)level=%d\n",level);
				 length+=8;
                 break;
			 }
			 break;
	}
	return (run<<24)|(level<<16)|(sign<<8)|length;
}




// draw the bitmap
void Display_Video(char* Text,HWND hWnd,int id,int XX,int YY,int W,int H,int Res,int Style,int Exstyle)
{
HWND A;
HBITMAP hBitmap;
MSG Msg;


// set default style
Style = WS_POPUPWINDOW|WS_VISIBLE|SS_BITMAP |WS_CAPTION;
Exstyle = WM_SETREDRAW|WM_PAINT;
// form for the image

A = CreateWindowEx(Exstyle,"Static","Decoded Video",Style,320,240,320,240,hWnd,NULL,NULL,NULL);


//ShowWindow(A,1);

//SetWindowPos(A,HWND_TOP,320,240,320,240,SWP_NOREDRAW);

// diplay decoded bmp sequence
do{


for(int i=start_frame-1;i<end_frame;i++){

hBitmap = CreateBitmap(320,240,1,32,i_buffer[i]);
SendMessage(A,(UINT)STM_SETIMAGE,(WPARAM)IMAGE_BITMAP,(LPARAM)hBitmap);
Sleep(33/speed_factor);
	
}

}while(is_loop);

DestroyWindow(A);


}


HWND WindowHandle(void)
{
HWND hConWnd;
OSVERSIONINFO os;
char szTempTitle[64], szClassName[128], szOriginalTitle[1024];
os.dwOSVersionInfoSize = sizeof( OSVERSIONINFO );
GetVersionEx( &os );
// may not work on WIN9x
if ( os.dwPlatformId == VER_PLATFORM_WIN32s ) return 0;
GetConsoleTitle( szOriginalTitle, sizeof ( szOriginalTitle ) );
sprintf( szTempTitle,"%u - %u", GetTickCount(), GetCurrentProcessId() );
SetConsoleTitle( szTempTitle );
Sleep( 40 );
// handle for NT and XP
hConWnd = FindWindow( NULL, szTempTitle );
SetConsoleTitle( szOriginalTitle );
// may not work on WIN9x
if ( os.dwPlatformId == VER_PLATFORM_WIN32_WINDOWS )
{
hConWnd = GetWindow( hConWnd, GW_CHILD );
if ( hConWnd == NULL ) return 0;
GetClassName( hConWnd, szClassName, sizeof ( szClassName ) );
// while ( _stricmp( szClassName, "ttyGrab" ) != 0 )
while ( strcmp( szClassName, "ttyGrab" ) != 0 )
{
hConWnd = GetNextWindow( hConWnd, GW_HWNDNEXT );
if ( hConWnd == NULL ) return 0;
GetClassName( hConWnd, szClassName, sizeof( szClassName ) );
}
}
return hConWnd;
}


int main(int argc,char **argv){

	
	printf("%15s","Start frame: ");
	scanf("%d",&start_frame);

	printf("%15s","End frame: ");
	scanf("%d",&end_frame);

	printf("%15s","Speed factor: ");
	scanf("%f",&speed_factor);

	printf("%15s","Loop(1/0)? ");
	scanf("%d",&is_loop);


	strcpy(outname,argv[1]);
	outname[strlen(outname)-4]=0;
	f = fopen(argv[1],"rb");

	fseek(f,0,SEEK_END);
	mpeg_size=ftell(f);
	rewind(f);
	
	



	mpgBuffer = (unsigned char*)malloc(sizeof(unsigned char)*mpeg_size);
	if (mpgBuffer==NULL)error("file error:out of memory!");
	fread(mpgBuffer,1,mpeg_size,f);
	rewind(f);
	fclose(f);

	hConWnd = WindowHandle();

	printf("\nDecode Start... ");
	decodeSequence();
   

	Display_Video(NULL,hConWnd,1,320,0,0,0);

	

	system("pause");
	return 0;
}
