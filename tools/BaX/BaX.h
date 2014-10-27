//////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                          //
//    BaX  Copyright (C) 2010-2012  Binyamin Chai                                           //
//    Published by "ibnatfatil"                                                             //
//    See the file AUTHORS for a complete list of people who contributed to the code        //
//                                                                                          //
//    This is BaX.h of the BaX-lib project.                                                 //
//    This is the only file you will need to add to your program in order to use BaX-lib.   //
//                                                                                          //
//    BaX is free software: you can redistribute it and/or modify                           //
//    it under the terms of the GNU General Public License as published by                  //
//    the Free Software Foundation, either version 3 of the License, or                     //
//    (at your option) any later version.                                                   //
//                                                                                          //
//    BaX is distributed in the hope that it will be useful,                                //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of                        //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                         //
//    GNU General Public License for more details.                                          //
//                                                                                          //
//    You should have received a copy of the GNU General Public License                     //
//    along with BaX.  If not, see <http://www.gnu.org/licenses/>.                          //
//                                                                                          //
//    BaX-lib version 1.1.2 pre                                                             //
//                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _libBaX_h_included_
#ifndef _BaX_h_included_
#define _BaX_h_included_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

using namespace std;


struct color
{
    unsigned char r,g,b;
    color()
    {
        r=g=b=0;
    }
    color(unsigned int hexval)
    {
        set(hexval);
    }
    color(unsigned char red,unsigned char green,unsigned char blue)
    {
        set(red,green,blue);
    }
    void set(color pix)
    {
        r=pix.r,g=pix.g,b=pix.b;
    }
    void setp(color pix)
    {
        r=pix.r%256,g=pix.g%256,b=pix.b%256;
    }
    void set(color *pix)
    {
        r=pix->r,g=pix->g,b=pix->b;
    }
    void set(unsigned char red,unsigned char green,unsigned char blue)
    {
        r=red,g=green,b=blue;
    }
    void set(unsigned int val)
    {
        r=(val/65536)%256,g=(val/256)%256,b=(val)%256;
    }
    void sethsv(long double hue,long double sat,long double val)
    {
        set(hsv(hue,sat,val));
    }
    void setscale(long double x,color pix1, color pix2)
    {
        set(scale(x,pix1,pix2));
    }
    unsigned int hsv(long double hue,long double sat,long double val)
    {
        hue/=60.0;
        val*=256.0;
        if(val>255.0) val=255.0;
        long double f=hue-floor(hue),p=val*(1.0-sat),q=val*(1.0-sat*f),t=val*(1.0-(1.0-f)*sat);
        switch((int)hue)
        {
            case 0: return chex(val,t,p); break;
            case 1: return chex(q,val,p); break;
            case 2: return chex(p,val,t); break;
            case 3: return chex(p,q,val); break;
            case 4: return chex(t,p,val); break;
            case 5: return chex(val,p,q); break;
            default: return 0; break;
        }
    }
    unsigned int scale(long double x,color pix1,color pix2)
    {
        unsigned char red=(unsigned char)abs((pix2.r-pix1.r)*x+pix1.r)%256;
        unsigned char green=(unsigned char)abs((pix2.g-pix1.g)*x+pix1.g)%256;
        unsigned char blue=(unsigned char)abs((pix2.b-pix1.b)*x+pix1.b)%256;
        return chex(red,green,blue);
    }
    bool operator==(color pix)
    {
        if(pix.r==r and pix.g==g and pix.b==b) return true;
        return false;
    }
    bool operator!=(color pix)
    {
        if(pix.r==r and pix.g==g and pix.b==b) return false;
        return true;
    }
    color operator!()
    {
        return (color)(0xFFFFFF-this->chex());
    }
    color operator/(color pix)
    {
        this->r=abs(pix.r-this->r);
        this->b=abs(pix.b-this->b);
        this->g=abs(pix.g-this->g);
        return *this;
    }
    unsigned int chex()
    {
        return 65536*r+256*g+b;
    }
    unsigned int chex(unsigned char red,unsigned char green,unsigned char blue)
    {
        return 65536*red+256*green+blue;
    }
};

class bitmap
{
protected:
    char fsize[4],width[4],height[4],rawsize[4],pfname[1024];
    unsigned char *raw;
    long int wsize,hsize,pads,row;
    struct floodstack
    {
	    struct list
	    {
		    long int x,y;
		    list *next;
	    };
	    list *start;
	    floodstack()
	    {
		    start=NULL;
	    }
	    ~floodstack()
	    {
		    list *temp;
		    while(start)
		    {
			    temp=start->next;
			    delete start;
			    start=temp;
		    }
	    }
	    bool set(long int x,long int y)
	    {
		    list *temp=start;
		    start=new(nothrow) list;
		    if(!start) return false;
		    start->x=x,start->y=y;
		    start->next=temp;
		    return true;
	    }
	    bool get(long int &x,long int &y)
	    {
		    if(start)
		    {
			    list *temp=start->next;
			    x=start->x,y=start->y;
			    delete start;
			    start=temp;
			    return true;
		    }
		    else return false;
	    }
    };
    long int padding(long int x)
    {
        const long int p=(3*x)%4;
        return p?(4-p):p;
    }
    bool dfslave(long int x,long int y,color *pix,long int max)
    {
        bool p=true;
        for(long int i=0;i<max;i++) if(pix[i]==get(x,y)) return false;
        return p;
    }
    bool colchk(long int x,long int y,bitmap *bmp,color *pix,long int n)
    {
        for(long int i=0;i<n;i++) if(bmp->upget(x,y)==pix[i]) return false;
        return true;
    }
public:
    static constexpr long double pi=3.14159265358979;
    unsigned long int foutd;
    bitmap()
    {
        cout<<"BaX-lib-error: empty bitmap cannot constructed";
        exit(1);
    }
    bitmap(bitmap &bmp)
    {
        foutd=bmp.foutd,wsize=bmp.wsize,hsize=bmp.hsize,pads=bmp.pads,row=bmp.row;
        memcpy(fsize,bmp.fsize,4),memcpy(width,bmp.width,4),memcpy(height,bmp.height,4),memcpy(rawsize,bmp.rawsize,4);
        raw=new(nothrow) unsigned char[le32touit(rawsize)];
        if(!raw)
        {
            cout<<"BaX-lib-error: allocating raw bitmap memory";
            exit(1);
        }
        memcpy(raw,bmp.raw,le32touit(rawsize));
    }
    bitmap(long int ws,long int hs)
    {
        foutd=0,wsize=ws,hsize=hs,pads=padding(wsize),row=3*wsize+pads;
        if(hs<1 or ws<1)
        {
            cout<<"BaX-lib-error: be reasonable demand the impossible";
            exit(1);
        }
        raw=new(nothrow) unsigned char[3*hsize*wsize+hsize*pads];
        if(!raw)
        {
            cout<<"BaX-lib-error: allocating raw bitmap memory";
            exit(1);
        }
        uitole32(fsize,3*hsize*wsize+hsize*pads+54);
        uitole32(width,wsize);
        uitole32(height,hsize);
        uitole32(rawsize,3*hsize*wsize+hsize*pads);
    }
    bitmap(long int ws,long int hs,color pix)
    {
        foutd=0,wsize=ws,hsize=hs,pads=padding(wsize),row=3*wsize+pads;
        if(hs<1 or ws<1)
        {
            cout<<"BaX-lib-error: be reasonable demand the impossible";
            exit(1);
        }
        raw=new(nothrow) unsigned char[3*hsize*wsize+hsize*pads];
        if(!raw)
        {
            cout<<"BaX-lib-error: allocating raw bitmap memory";
            exit(1);
        }
        uitole32(fsize,3*hsize*wsize+hsize*pads+54);
        uitole32(width,wsize);
        uitole32(height,hsize);
        uitole32(rawsize,3*hsize*wsize+hsize*pads);
        fillup(pix);
    }
    bitmap(char *filename)
    {
        char magic[2],dummy[4];
		ifstream fin(filename,ios::binary);
		if(!fin)
		{
			cout<<"BaX-lib-error: unable to open \""<<filename<<"\""<<endl;
			exit(1);
		}
        fin.read(magic,2);
        if(!(magic[0]=='B' and magic[1]=='M')) cout<<"BaX-lib-warning: input file seems unsupported because the magicnumber is \""<<magic<<"\""<<endl;
        fin.seekg(0x1C,ios::beg);
        fin.read(magic,2);
        if(!(magic[0]==24 and magic[1]==0))
        {
            cout<<"BaX-lib-error: only 24bit input files are supported and this is a "<<le32touit(magic)<<"bit bitmap"<<endl;
            fin.close();
            exit(1);
        }
        fin.seekg(0x12,ios::beg);
        fin.read(width,4);
        fin.seekg(0x16,ios::beg);
        fin.read(height,4);
        foutd=0,wsize=le32touit(width),hsize=le32touit(height),pads=padding(wsize),row=3*wsize+pads;
        if(wsize<1 or hsize<1)
        {
            cout<<"BaX-lib-error: unreal 24bit bitmap file"<<endl;
            fin.close();
            exit(1);
        }
        fin.seekg(0x0A,ios::beg);
        fin.read(dummy,4);
        int start=le32touit(dummy);
        fin.seekg(start,ios::beg);
        raw=new (nothrow) unsigned char[3*hsize*wsize+hsize*pads];
        if(!raw)
        {
            cout<<"BaX-lib-error: allocating raw bitmap memory"<<endl;
            fin.close();
            exit(1);
        }
        fin.read((char*)raw,3*hsize*wsize+hsize*pads);
        fin.close();
        uitole32(fsize,3*hsize*wsize+hsize*pads+54);
        uitole32(rawsize,3*hsize*wsize+hsize*pads);
    }
    ~bitmap()
    {
        delete[] raw;
    }
    bool prtct(long int x,long int y)
    {
        return (x>=0 and x<wsize and y>=0 and y<hsize);
    }
    void protect(long int &x,long int &y)
    {
        if(x>=wsize) x%=wsize;
        else if(x<0 and abs(x)<wsize) x=wsize-abs(x)-1;
        if(y>=hsize) y%=hsize;
        else if(y<0 and abs(y)<hsize) y=hsize-abs(y)-1;
    }
    color upget(long int x,long int y)
    {
        color pix;
        const long int a=row*(hsize-y-1)+3*x;
        pix.b=raw[a],pix.g=raw[a+1],pix.r=raw[a+2];
        return pix;
    }
    bool upset(long int x,long int y,color pix)
    {
        const long int a=row*(hsize-y-1)+3*x;
        raw[a]=pix.b,raw[a+1]=pix.g,raw[a+2]=pix.r;
        return true;
    }
    color get(long int x,long int y)
    {
        color pix;
        protect(x,y);
        const long int a=row*(hsize-y-1)+3*x;
        pix.b=raw[a],pix.g=raw[a+1],pix.r=raw[a+2];
        return pix;
    }
    bool set(long int x,long int y,color pix)
    {
        protect(x,y);
        const long int a=row*(hsize-y-1)+3*x;
        raw[a]=pix.b,raw[a+1]=pix.g,raw[a+2]=pix.r;
        return true;
    }
    void set(long int fat,long int x, long int y,color pix)
    {
        for(long int i=0;i<fat;i++) for(long int j=0;j<fat;j++)
        {
            set(x+i,y+j,pix);
            set(x-i,y-j,pix);
            set(x+i,y-j,pix);
            set(x-i,y+j,pix);
        }
        set(x,y,pix);
    }
    void set(long int fat,long int x, long int y,color pix,color pixsh)
    {
        for(long int i=0;i<fat;i++) for(long int j=0;j<fat;j++)
        {
            if(pix!=get(x+i,y+j)) set(x+i,y+j,pixsh);
            if(pix!=get(x-i,y-j)) set(x-i,y-j,pixsh);
            if(pix!=get(x+i,y-j)) set(x+i,y-j,pixsh);
            if(pix!=get(x-i,y+j)) set(x-i,y+j,pixsh);
        }
        set(x,y,pix);
    }
    color subget(long int binx,long int biny,long int x,long int y)
    {
        long int sum=0;
        for(long int i=0;i<binx;i++) for(long int j=0;j<biny;j++) sum+=gcol(x+i,y+j);
        return (color)(sum/(binx*biny));
    }
    void subset(long int binx,long int biny,long int x,long int y,color pix)
    {
        for(long int i=0;i<binx;i++) for(long int j=0;j<biny;j++) set(x+i,y+j,pix);
    }
    long int gcol(long int x,long int y)
    {
        return get(x,y).chex();
    }
    long int upgcol(long int x,long int y)
    {
        return upget(x,y).chex();
    }
    void chcol(color in,color out)
    {
        for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++) if(in==get(i,j)) upset(i,j,out);
    }
    void fillup(color pix)
    {
        for(long int x=0;x<wsize;x++) for(long int y=0;y<hsize;y++) upset(x,y,pix);
    }
	void inline csnow(color *val,long int max)
	{
		for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++) upset(i,j,val[rand()%max]);
	}
	void snow(color *val,long int max)
	{
        bool p;
		for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++)
		{
            p=true;
			for(long int k=0;k<max;k++) if(upget(i,j)==val[k])
            {
                p=false;
                break;
            }
			if(p) upset(i,j,rand()%16777216);
		}
	}
    void flipx()
    {
        for(long int j=0;j<hsize;j++) for(long int i=0;i<(wsize/2);i++) swap(i,j,wsize-i-1,j);
    }
    void flipy()
    {
        for(long int i=0;i<wsize;i++) for(long int j=0;j<(hsize/2);j++) swap(i,j,i,hsize-j-1);
    }
    void rotaclock()
    {
        if(wsize!=hsize) return;
        bitmap *bmp=new(nothrow) bitmap(wsize,hsize);
        bmp->import(this);
        for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++) upset(i,j,bmp->get(-j,i));
        delete bmp;
    }
    void rotclock()
    {
        if(wsize!=hsize) return;
        bitmap *bmp=new(nothrow) bitmap(wsize,hsize);
        bmp->import(this);
        for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++) upset(i,j,bmp->get(j,-i));
        delete bmp;
    }
    void rotate(long double t)
    {
        if(wsize!=hsize) return;
        bitmap *bmp=new(nothrow) bitmap(wsize,hsize);
        bmp->import(this);
        for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++)
        upset(i,j,bmp->get((long int)(cos(t)*i-sin(t)*j+wsize)%wsize,(long int)(cos(t)*i+sin(t)*j+hsize)%hsize));
        delete bmp;
    }
    void negative()
    {
        for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++) upset(i,j,!get(i,j));
    }
    void negative(color pix)
    {
        for(long int i=0;i<wsize;i++) for(long int j=0;j<hsize;j++) upset(i,j,get(i,j)/pix);
    }
    void swap(long int x1,long int y1,long int x2,long int y2)
    {
        color pix=get(x1,y1);
        set(x1,y1,get(x2,y2));
        set(x2,y2,pix);
    }
    void param(long int (*x)(long double),long int (*y)(long double),color (*pix)(long double))
    {
        const long int dx=x(1.0)-x(0.0),dy=y(1.0)-y(0.0);
        const long double dt=(abs(dx)>abs(dy))?(1/(long double)(abs(dx)+1)):(1/(long double)(abs(dy)+1));
        for(long double t=0.0;t<=1.0;t+=dt) set(x(t),y(t),pix(t));
    }
    void param(long int fat,long int (*x)(long double),long int (*y)(long double),color (*pix)(long double))
    {
        const long int dx=x(1.0)-x(0.0),dy=y(1.0)-y(0.0);
        const long double dt=(abs(dx)>abs(dy))?(1/(long double)(abs(dx)+1)):(1/(long double)(abs(dy)+1));
        for(long double t=0.0;t<=1.0;t+=dt) set(fat,x(t),y(t),pix(t));
    }
    void param(long int fat,long int (*x)(long double),long int (*y)(long double),color (*pix)(long double),color (*pixsh)(long double))
    {
        const long int dx=x(1.0)-x(0.0),dy=y(1.0)-y(0.0);
        const long double dt=(abs(dx)>abs(dy))?(1/(long double)(abs(dx)+1)):(1/(long double)(abs(dy)+1));
        for(long double t=0.0;t<=1.0;t+=dt) set(fat,x(t),y(t),pix(t),pixsh(t));
    }
    void connex(long int x1,long int y1,long int x2,long int y2,color pix)
    {
        const long int dx=x2-x1,dy=y2-y1;
        if(!(dx==0 and dy==0))
        {
            const long double dt=(abs(dx)>abs(dy))?(1/(long double)(abs(dx)+1)):(1/(long double)(abs(dy)+1));
            for(long double t=0.0;t<=1.0;t+=dt) set(x1+dx*t,y1+dy*t,pix);
        }
        else set(x1,y1,pix);
    }
    void connex(long int fat,long int x1,long int y1,long int x2,long int y2,color pix)
    {
        const long int dx=x2-x1,dy=y2-y1;
        if(!(dx==0 and dy==0))
        {
            const long double dt=(abs(dx)>abs(dy))?(1/(long double)(abs(dx)+1)):(1/(long double)(abs(dy)+1));
            for(long double t=0.0;t<=1.0;t+=dt) set(fat,x1+dx*t,y1+dy*t,pix);
        }
        else set(fat,x1,y1,pix);
    }
    void connex(long int fat,long int x1,long int y1,long int x2,long int y2,color pix,color pixsh)
    {
        const long int dx=x2-x1,dy=y2-y1;
        if(!(dx==0 and dy==0))
        {
            const long double dt=(abs(dx)>abs(dy))?(1/(long double)(abs(dx)+1)):(1/(long double)(abs(dy)+1));
            for(long double t=0.0;t<=1.0;t+=dt) set(fat,x1+dx*t,y1+dy*t,pix,pixsh);
        }
        else set(fat,x1,y1,pix,pixsh);
    }
	void box(long int x1,long int y1,long int x2,long int y2,color pix)
	{
		connex(x1,y1,x2,y1,pix);
		connex(x1,y1,x1,y2,pix);
		connex(x2,y1,x2,y2,pix);
		connex(x1,y2,x2,y2,pix);
	}
	void box(long int fat,long int x1,long int y1,long int x2,long int y2,color pix)
	{
		connex(fat,x1,y1,x2,y1,pix);
		connex(fat,x1,y1,x1,y2,pix);
		connex(fat,x2,y1,x2,y2,pix);
		connex(fat,x1,y2,x2,y2,pix);
	}
	void box(long int fat, long int x1,long int y1,long int x2,long int y2,color pix,color pixsh)
	{
		connex(fat,x1,y1,x2,y1,pix,pixsh);
		connex(fat,x1,y1,x1,y2,pix,pixsh);
		connex(fat,x2,y1,x2,y2,pix,pixsh);
		connex(fat,x1,y2,x2,y2,pix,pixsh);
	}
	void solidbox(long int x1,long int y1,long int x2,long int y2,color pix)
	{
	    subset(abs(x2-x1),abs(y2-y1),x1<x2?x1:x2,y1<y2?y1:y2,pix);
	}
	void solidbox(long int x1,long int y1,long int x2,long int y2,color pix,color pixsh)
	{
	    solidbox(x1,y1,x2,y2,pix);
	    box(x1,y1,x2,y2,pixsh);
	}
    void solidbox(long int fat, long int x1,long int y1,long int x2,long int y2,color pix,color pixsh)
	{
	    solidbox(x1,y1,x2,y2,pix);
	    box(fat,x1+(x1<x2?fat-1:-1*fat+1),y1+(y1<y2?fat-1:-1*fat+1),x2-(x1<x2?fat-1:-1*fat+1),y2-+(y1<y2?fat-1:-1*fat+1),pixsh);
	}
	void vector(long int x,long int y,long double r,long double theta,color pix)
	{
		connex(x,y,(r*cos(theta)+x),(r*sin(theta)+y),pix);
	}
	void vector(long int fat,long int x,long int y,long double r,long double theta,color pix)
	{
		connex(fat,x,y,(r*cos(theta)+x),(r*sin(theta)+y),pix);
	}
	void vector(long int fat,long int x,long int y,long double r,long double theta,color pix,color pixsh)
	{
		connex(fat,x,y,(r*cos(theta)+x),(r*sin(theta)+y),pix,pixsh);
	}
	void kykel(long double r,long int x,long int y,color pix)
	{
		for(long double t=0.0;t<2.0*pi;t+=1.0/(2.0*pi*r+1.0)) set(round(r*sin(t)+x),round(r*cos(t)+y),pix);
	}
	void kykel(long int fat,long double r,long int x,long int y,color pix)
	{
		for(long int i=-1*fat;i<=fat;i++) for(long double t=0.0;t<2.0*pi;t+=1.0/(2.0*pi*r+1.0)) set(round((r+i)*sin(t)+x),round((r+i)*cos(t)+y),pix);
	}
	void kykel(long int fat,long double r,long int x,long int y,color pix,color pixsh)
	{
         kykel(fat,r,x,y,pixsh);
         kykel(r,x,y,pix);
	}
    void solidkykel(long double r,long int x,long int y,color pix)
    {
        for(long int i=0;i<=r;i++) kykel(i,x,y,pix);
    }
    void solidkykel(long double r,long int x,long int y,color pix,color pixsh)
    {
        solidkykel(r-1,x,y,pix);
        kykel(r,x,y,pixsh);
    }
    void solidkykel(long int fat,long double r,long int x,long int y,color pix,color pixsh)
    {
        solidkykel(abs(r-fat),x,y,pix);
        for(long int i=0;i<=fat;i++) kykel(abs(r-fat)+i,x,y,pixsh);
    }
	void elkyk(long double rx,long double ry,long int x,long int y,color pix)
	{
		for(long double t=0.0;t<2.0*pi;t+=(rx>ry)?1.0/(2.0*pi*rx+1.0):1.0/(2.0*pi*ry+1.0)) set(rx*sin(t)+x,ry*cos(t)+y,pix);
	}
	void elkyk(long int fat,long double rx,long double ry,long int x,long int y,color pix)
	{
		for(long int i=-1*fat;i<=fat;i++) for(long double t=0.0;t<2.0*pi;t+=(rx>ry)?1.0/(2.0*pi*rx+1.0):1.0/(2.0*pi*ry+1.0))
        set((rx+i)*sin(t)+x,(ry+i)*cos(t)+y,pix);
	}
	void elkyk(long int fat,long double rx,long double ry,long int x,long int y,color pix,color pixsh)
	{
         elkyk(fat,rx,ry,x,y,pixsh);
         elkyk(rx,ry,x,y,pix);
	}
	void solidelkyk(long double rx,long double ry,long int x,long int y,color pix)
	{
	    for(long int i=(rx<ry?rx:ry);i>=0;i--) elkyk(rx-i,ry-i,x,y,pix);
	}
	void solidelkyk(long double rx,long double ry,long int x,long int y,color pix,color pixsh)
	{
        solidelkyk(rx,ry,x,y,pix);
        elkyk(rx,ry,x,y,pixsh);
	}
	void solidelkyk(long int fat, long double rx,long double ry,long int x,long int y,color pix,color pixsh)
	{
        solidelkyk(rx,ry,x,y,pix);
        elkyk(fat,rx-fat,ry-fat,x,y,pixsh);
	}
    void pc(long int *x,long int *y,long int length,color pix)
    {
        for(long int i=0;i<length-1;i++) connex(x[i],y[i],x[i+1],y[i+1],pix);
    }
    void pc(long int fat,long int *x,long int *y,long int length,color pix)
    {
        for(long int i=0;i<length-1;i++) connex(fat,x[i],y[i],x[i+1],y[i+1],pix);
    }
    void pc(long int fat,long int *x,long int *y,long int length,color pix,color pixsh)
    {
        for(long int i=0;i<length-1;i++) connex(fat,x[i],y[i],x[i+1],y[i+1],pix,pixsh);
    }
    void flood4(long int x,long int y,color pix)
	{
        color delim=get(x,y);
		floodstack execute;
		while(true)
		{
			set(x,y,pix);
			if(delim==get(x+1,y)) execute.set(x+1,y);
			if(delim==get(x-1,y)) execute.set(x-1,y);
			if(delim==get(x,y+1)) execute.set(x,y+1);
			if(delim==get(x,y-1)) execute.set(x,y-1);
			if(!execute.get(x,y)) break;
		}
	}
	void flood8(long int x,long int y,color pix)
	{
        color delim=get(x,y);
		floodstack execute;
		while(true)
		{
			set(x,y,pix);
			if(delim==get(x+1,y)) execute.set(x+1,y);
			if(delim==get(x-1,y)) execute.set(x-1,y);
			if(delim==get(x,y+1)) execute.set(x,y+1);
			if(delim==get(x,y-1)) execute.set(x,y-1);
			if(delim==get(x+1,y+1)) execute.set(x+1,y+1);
			if(delim==get(x-1,y-1)) execute.set(x-1,y-1);
			if(delim==get(x-1,y+1)) execute.set(x-1,y+1);
			if(delim==get(x+1,y-1)) execute.set(x+1,y-1);
			if(!execute.get(x,y)) break;
		}
	}
    void dflood4(long int x,long int y,color pix)
	{
		floodstack execute;
		while(true)
		{
			set(x,y,pix);
			if(pix!=get(x+1,y) and prtct(x+1,y)) execute.set(x+1,y);
			if(pix!=get(x-1,y) and prtct(x-1,y)) execute.set(x-1,y);
			if(pix!=get(x,y+1) and prtct(x,y+1)) execute.set(x,y+1);
			if(pix!=get(x,y-1) and prtct(x,y-1)) execute.set(x,y-1);
			if(!execute.get(x,y)) break;
		}
	}
	void dflood8(long int x,long int y,color pix)
	{
		floodstack execute;
		while(true)
		{
			set(x,y,pix);
			if(pix!=get(x+1,y) and prtct(x+1,y)) execute.set(x+1,y);
			if(pix!=get(x-1,y) and prtct(x-1,y)) execute.set(x-1,y);
			if(pix!=get(x,y+1) and prtct(x,y+1)) execute.set(x,y+1);
			if(pix!=get(x,y-1) and prtct(x,y-1)) execute.set(x,y-1);
			if(pix!=get(x+1,y+1) and prtct(x+1,y+1)) execute.set(x+1,y+1);
			if(pix!=get(x-1,y-1) and prtct(x-1,y-1)) execute.set(x-1,y-1);
			if(pix!=get(x-1,y+1) and prtct(x-1,y+1)) execute.set(x-1,y+1);
			if(pix!=get(x+1,y-1) and prtct(x+1,y-1)) execute.set(x+1,y-1);
			if(!execute.get(x,y)) break;
		}
	}
    void dflood4(long int x,long int y,color *pix,long int n)
	{
		floodstack execute;
		while(true)
		{
			set(x,y,pix[0]);
			if(dfslave(x+1,y,pix,n) and prtct(x+1,y)) execute.set(x+1,y);
			if(dfslave(x-1,y,pix,n) and prtct(x-1,y)) execute.set(x-1,y);
			if(dfslave(x,y+1,pix,n) and prtct(x,y+1)) execute.set(x,y+1);
			if(dfslave(x,y-1,pix,n) and prtct(x,y-1)) execute.set(x,y-1);
			if(!execute.get(x,y)) break;
		}
	}
	void dflood8(long int x,long int y,color *pix,long int n)
	{
		floodstack execute;
		while(true)
		{
			set(x,y,pix[0]);
			if(dfslave(x+1,y,pix,n) and prtct(x+1,y)) execute.set(x+1,y);
			if(dfslave(x-1,y,pix,n) and prtct(x-1,y)) execute.set(x-1,y);
			if(dfslave(x,y+1,pix,n) and prtct(x,y+1)) execute.set(x,y+1);
			if(dfslave(x,y-1,pix,n) and prtct(x,y-1)) execute.set(x,y-1);
			if(dfslave(x+1,y+1,pix,n) and prtct(x+1,y+1)) execute.set(x+1,y+1);
			if(dfslave(x-1,y-1,pix,n) and prtct(x-1,y-1)) execute.set(x-1,y-1);
			if(dfslave(x-1,y+1,pix,n) and prtct(x-1,y+1)) execute.set(x-1,y+1);
			if(dfslave(x+1,y-1,pix,n) and prtct(x+1,y-1)) execute.set(x+1,y-1);
			if(!execute.get(x,y)) break;
		}
	}
    bool compare(bitmap *bmp)
    {
        if(bmp->cwsize()==wsize and bmp->chsize()==hsize and bmp->crawsize()==crawsize()) return (memcmp(bmp->raw,raw,crawsize()))==0?true:false;
        else return false;
    }
    bool swap(bitmap *bmp)
    {
        if(bmp->cwsize()==wsize and bmp->chsize()==hsize and bmp->crawsize()==crawsize())
        {
            unsigned char *temp=bmp->raw;
            bmp->raw=raw;
            raw=temp;
            return true;
        }
        else return false;
    }
    bool import(bitmap *bmp)
    {
        if(bmp->crawsize()==crawsize() and bmp->cwsize()==wsize and bmp->chsize()==hsize)
        {
            memcpy(raw,bmp->raw,crawsize());
            return true;
        }
        else return false;
    }
    void import(long int x, long int y,bitmap *bmp)
    {
        for(long int i=0;i<bmp->cwsize();i++) for(long int j=0;j<bmp->chsize();j++) set(x+i,y+j,bmp->upget(i,j));
    }
    void import(long int x, long int y,bitmap *bmp,color *pix,long int n)
    {
        for(long int i=0;i<bmp->cwsize();i++) for(long int j=0;j<bmp->chsize();j++) if(colchk(i,j,bmp,pix,n)) set(x+i,y+j,bmp->upget(i,j));
    }
    bool save(char *fname)
    {
        if(pfname!=fname) strcpy(pfname,fname);
        ofstream fout(fname,ios::binary);
        if(!fout) return false;
        static const char magicnr[2]={66,77},aoh[12]={66,67,118,90,54,0,0,0,40,0,0,0};
        static const char nbt[8]={1,0,24,0,0,0,0,0},whnn[16]={-100,14,0,0,-100,14,0,0,0,0,0,0,0,0,0,0};
        fout.write(magicnr,2);
        fout.write(fsize,4);
        fout.write(aoh,12);
        fout.write(width,4);
        fout.write(height,4);
        fout.write(nbt,8);
        fout.write(rawsize,4);
        fout.write(whnn,16);
        fout.write((char*)raw,3*hsize*wsize+hsize*pads);
        fout.close();
        foutd++;
        return true;
    }
    bool save_wx(char *fname)
    {
        if(fname!=pfname) strcpy(pfname,fname);
        strcat(pfname,".bmp");
        return save(pfname);
    }
    bool save_wx_seq(char *fname)
    {
        if(fname!=pfname) strcpy(pfname,fname);
        strcat(pfname,itoa(foutd));
        return save_wx(pfname);
    }
    bool save_jpeg(char *fname)
    {
        if(!system("mogrify -version > /dev/null"))
        {
            bool p,q;
            static char str[1024];
            p=save_wx(fname);
            strcpy(str,"mogrify -format jpeg ");
            strcat(str,pfname);
            strcat(str," && rm ");
            strcat(str,pfname);
            strcat(str," &");
	        q=system(str);
	        if(p and q) return true;
	        else return false;
	    }
	    else return false;
    }
    bool save_jpeg_seq(char *fname)
    {
        if(!system("mogrify -version > /dev/null")) //some systems require a ! here
        {
            bool p,q;
            static char str[1024];
            p=save_wx_seq(fname);
            strcpy(str,"mogrify -format jpeg ");
            strcat(str,pfname);
            strcat(str," && rm ");
            strcat(str,pfname);
            strcat(str," &");
	        q=system(str);
	        if(p and q) return true;
	        else return false;
	    }
	    else return false;
    }
    bool save_png(char *fname)
    {
        if(!system("mogrify -version > /dev/null"))
        {
            bool p,q;
            static char str[1024];
            p=save_wx(fname);
            strcpy(str,"mogrify -format png ");
            strcat(str,pfname);
            strcat(str," && rm ");
            strcat(str,pfname);
            strcat(str," &");
	        q=system(str);
	        if(p and q) return true;
	        else return false;
	    }
	    else return false;
    }
    bool save_png_seq(char *fname)
    {
        if(!system("mogrify -version > /dev/null")) //some systems require a ! here
        {
            bool p,q;
            static char str[1024];
            p=save_wx_seq(fname);
            strcpy(str,"mogrify -format png ");
            strcat(str,pfname);
            strcat(str," && rm ");
            strcat(str,pfname);
            strcat(str," &");
	        q=system(str);
	        if(p and q) return true;
	        else return false;
	    }
	    else return false;
    }
    char *fname()
    {
        return foutd?pfname:NULL;
    }
    void uitole32(char le32[4],long int a)
    {
        le32[0]=(a%256),le32[1]=(a/256)%256,le32[2]=(a/65536)%256,le32[3]=(a/16777216)%256;
	}
    long int le32touit(char le32[4])
    {
        return (unsigned char)le32[0]+256*(unsigned char)le32[1]+65536*(unsigned char)le32[2]+16777216*(unsigned char)le32[3];
    }
    char *itoa(long int i)
    {
        static const int size=(8*sizeof(long int)*0.301029996)+2;
        static char a[size];
        sprintf(a,"%ld",i);
        return a;
    }
    void randinit()
    {
        srand(time(NULL));
    }
    long int cwsize()
    {
        return wsize;
    }
    long int chsize()
    {
        return hsize;
    }
    long int csize()
    {
        return le32touit(fsize);
    }
    long int crawsize()
    {
        return le32touit(rawsize);
    }
    bool craw(unsigned char *buffer,long int size)
    {
        if(size==crawsize())
        {
            memcpy(buffer,raw,size);
            return true;
        }
        else return false;
    }
};

const color BLACK(0x000000);
const color NAVY(0x000080);
const color DARKBLUE(0x00008B);
const color MEDIUMBLUE(0x0000CD);
const color BLUE(0x0000FF);
const color DARKGREEN(0x006400);
const color GREEN(0x008000);
const color TEAL(0x008080);
const color DARKCYAN(0x008B8B);
const color DEEPSKYBLUE(0x00BFFF);
const color DARKTURQUOISE(0x00CED1);
const color MEDIUMSPRINGGREEN(0x00FA9A);
const color LIME(0x00FF00);
const color SPRINGGREEN(0x00FF7F);
const color AQUA(0x00FFFF);
const color CYAN(0x00FFFF);
const color MIDNIGHTBLUE(0x191970);
const color DODGERBLUE(0x1E90FF);
const color LIGHTSEAGREEN(0x20B2AA);
const color FORESTGREEN(0x228B22);
const color SEAGREEN(0x2E8B57);
const color DARKSLATEGRAY(0x2F4F4F);
const color LIMEGREEN(0x32CD32);
const color MEDIUMSEAGREEN(0x3CB371);
const color TURQUOISE(0x40E0D0);
const color ROYALBLUE(0x4169E1);
const color STEELBLUE(0x4682B4);
const color DARKSLATEBLUE(0x483D8B);
const color MEDIUMTURQUOISE(0x48D1CC);
const color INDIGO(0x4B0082);
const color DARKOLIVEGREEN(0x556B2F);
const color CADETBLUE(0x5F9EA0);
const color CORNFLOWERBLUE(0x6495ED);
const color MEDIUMAQUAMARINE(0x66CDAA);
const color DIMGRAY(0x696969);
const color SLATEBLUE(0x6A5ACD);
const color OLIVEDRAB(0x6B8E23);
const color SLATEGRAY(0x708090);
const color LIGHTSLATEGRAY(0x778899);
const color MEDIUMSLATEBLUE(0x7B68EE);
const color LAWNGREEN(0x7CFC00);
const color CHARTREUSE(0x7FFF00);
const color AQUAMARINE(0x7FFFD4);
const color MAROON(0x800000);
const color PURPLE(0x800080);
const color OLIVE(0x808000);
const color GRAY(0x808080);
const color SKYBLUE(0x87CEEB);
const color LIGHTSKYBLUE(0x87CEFA);
const color BLUEVIOLET(0x8A2BE2);
const color DARKRED(0x8B0000);
const color DARKMAGENTA(0x8B008B);
const color SADDLEBROWN(0x8B4513);
const color DARKSEAGREEN(0x8FBC8F);
const color LIGHTGREEN(0x90EE90);
const color MEDIUMPURPLE(0x9370D8);
const color DARKVIOLET(0x9400D3);
const color PALEGREEN(0x98FB98);
const color DARKORCHID(0x9932CC);
const color YELLOWGREEN(0x9ACD32);
const color SIENNA(0xA0522D);
const color BROWN(0xA52A2A);
const color DARKGRAY(0xA9A9A9);
const color LIGHTBLUE(0xADD8E6);
const color GREENYELLOW(0xADFF2F);
const color PALETURQUOISE(0xAFEEEE);
const color LIGHTSTEELBLUE(0xB0C4DE);
const color POWDERBLUE(0xB0E0E6);
const color FIREBRICK(0xB22222);
const color DARKGOLDENROD(0xB8860B);
const color MEDIUMORCHID(0xBA55D3);
const color ROSYBROWN(0xBC8F8F);
const color DARKKHAKI(0xBDB76B);
const color SILVER(0xC0C0C0);
const color MEDIUMVIOLETRED(0xC71585);
const color INDIANRED(0xCD5C5C);
const color PERU(0xCD853F);
const color CHOCOLATE(0xD2691E);
const color TAN(0xD2B48C);
const color LIGHTGREY(0xD3D3D3);
const color PALEVIOLETRED(0xD87093);
const color THISTLE(0xD8BFD8);
const color ORCHID(0xDA70D6);
const color GOLDENROD(0xDAA520);
const color CRIMSON(0xDC143C);
const color GAINSBORO(0xDCDCDC);
const color PLUM(0xDDA0DD);
const color BURLYWOOD(0xDEB887);
const color LIGHTCYAN(0xE0FFFF);
const color LAVENDER(0xE6E6FA);
const color DARKSALMON(0xE9967A);
const color VIOLET(0xEE82EE);
const color PALEGOLDENROD(0xEEE8AA);
const color LIGHTCORAL(0xF08080);
const color KHAKI(0xF0E68C);
const color ALICEBLUE(0xF0F8FF);
const color HONEYDEW(0xF0FFF0);
const color AZURE(0xF0FFFF);
const color SANDYBROWN(0xF4A460);
const color WHEAT(0xF5DEB3);
const color BEIGE(0xF5F5DC);
const color WHITESMOKE(0xF5F5F5);
const color MINTCREAM(0xF5FFFA);
const color GHOSTWHITE(0xF8F8FF);
const color SALMON(0xFA8072);
const color ANTIQUEWHITE(0xFAEBD7);
const color LINEN(0xFAF0E6);
const color LIGHTGOLDENRODYELLOW(0xFAFAD2);
const color OLDLACE(0xFDF5E6);
const color RED(0xFF0000);
const color FUCHSIA(0xFF00FF);
const color MAGENTA(0xFF00FF);
const color DEEPPINK(0xFF1493);
const color ORANGERED(0xFF4500);
const color TOMATO(0xFF6347);
const color HOTPINK(0xFF69B4);
const color CORAL(0xFF7F50);
const color DARKORANGE(0xFF8C00);
const color LIGHTSALMON(0xFFA07A);
const color ORANGE(0xFFA500);
const color LIGHTPINK(0xFFB6C1);
const color PINK(0xFFC0CB);
const color GOLD(0xFFD700);
const color PEACHPUFF(0xFFDAB9);
const color NAVAJOWHITE(0xFFDEAD);
const color MOCCASIN(0xFFE4B5);
const color BISQUE(0xFFE4C4);
const color MISTYROSE(0xFFE4E1);
const color BLANCHEDALMOND(0xFFEBCD);
const color PAPAYAWHIP(0xFFEFD5);
const color LAVENDERBLUSH(0xFFF0F5);
const color SEASHELL(0xFFF5EE);
const color CORNSILK(0xFFF8DC);
const color LEMONCHIFFON(0xFFFACD);
const color FLORALWHITE(0xFFFAF0);
const color SNOW(0xFFFAFA);
const color YELLOW(0xFFFF00);
const color LIGHTYELLOW(0xFFFFE0);
const color IVORY(0xFFFFF0);
const color WHITE(0xFFFFFF);

#endif
#endif

