/********************************************************************

Copyright (C), 2019, All rights reserved

File Name     :    UVolume.h
Description   :
History       :

<author>            <time>            <desc>
Ang Li             2019/6/8           create

********************************************************************/
#ifndef __VOLUME_H
#define __VOLUME_H

#include <iostream>

//A template class for 3-dimension value
template <class T>
class UVolume3D
{
public:
	T x;
	T y;
	T z;
public:
	UVolume3D() { x = 0; y = 0; z = 0; };
	UVolume3D(T xx, T yy, T zz)
	{
		x = xx; y = yy; z = zz;
	}
	UVolume3D(const UVolume3D & temp)
	{
		x = temp.x;
		y = temp.y;
		z = temp.z;
	}
	void set(T xx, T yy, T zz)
	{
		x = xx;
		y = yy;
		z = zz;
	}
	void exchange(UVolume3D &src)//���������н���
	{
		UVolume3D temp = src;
		src = *this;
		*this = temp;
	}
	UVolume3D crossMul(UVolume3D &vector)const//������ˣ�����
	{
		UVolume3D result;
		result.x = y * vector.z - vector.y * z;
		result.y = z * vector.x - vector.z * x;
		result.z = x * vector.y - vector.x * y;
		return result;
	}
	//���������
	UVolume3D &operator=(const UVolume3D & temp)
	{
		if (this != &temp)
		{
			x = temp.x;
			y = temp.y;
			z = temp.z;
		}
		return *this;
	}
	UVolume3D operator+(const UVolume3D &temp)const// '+'
	{
		return UVolume3D(x + temp.x, y + temp.y, z + temp.z);
	}
	UVolume3D operator-(const UVolume3D &temp)const// '-'
	{
		return UVolume3D(x - temp.x, y - temp.y, z - temp.z);
	}
	UVolume3D operator+(double k)const// '+'//20190604 updated
	{
		return UVolume3D(x + k, y + k, z + k);
	}
	UVolume3D operator-(double k)const// '-'
	{
		return UVolume3D(x - k, y - k, z - k);
	}
	T operator*(const UVolume3D &temp)const// '*'���
	{
		return x * temp.x + y * temp.y + z * temp.z;
	}
	UVolume3D operator*(double k)const// '*k'
	{
		return UVolume3D(x*k, y*k, z*k);
	}
	UVolume3D operator/(double k)const// '/k'
	{
		return UVolume3D(x / k, y / k, z / k);
	}
	//������ģ���࣬��ȷ���������ͣ��޷�ʹ�ñ�׼�����ֻ����cout
	friend std::ostream &operator<<(std::ostream &output, const UVolume3D &temp)//friend��<<���������Ϊ��Ԫ����
	{

		output << "( " << temp.x << " , " << temp.y << " , " << temp.z << " ) ";
		return output;
	}
};

template <class T>
class UVolume2D
{
public:
	T x;
	T y;
public:
	UVolume2D() { x = 0; y = 0; };
	UVolume2D(T xx, T yy)
	{
		x = xx; y = yy;
	}
	UVolume2D(const UVolume2D & temp)
	{
		x = temp.x;
		y = temp.y;
	}
	void set(T xx, T yy)
	{
		x = xx;
		y = yy;
	}
	void exchange(UVolume2D &src)
	{
		UVolume2D temp = src;
		src = *this;
		*this = temp;
	}
	UVolume2D &operator=(const UVolume2D & temp)
	{
		if (this != &temp)
		{
			x = temp.x;
			y = temp.y;
		}
		return *this;
	}
	UVolume2D operator+(const UVolume2D &temp)const// '+'
	{
		return UVolume2D(x + temp.x, y + temp.y);
	}
	UVolume2D operator-(const UVolume2D &temp)const// '-'
	{
		return UVolume2D(x - temp.x, y - temp.y);
	}
	UVolume2D operator+(double k)const// '+'//20190604 updated
	{
		return UVolume2D(x + k, y + k);
	}
	UVolume2D operator-(double k)const// '-'
	{
		return UVolume2D(x - k, y - k);
	}
	T operator*(const UVolume2D &temp)const// '*'���
	{
		return x * temp.x + y * temp.y;
	}
	UVolume2D operator*(double k)const// '*k'
	{
		return UVolume2D(x*k, y*k);
	}
	UVolume2D operator/(double k)const// '/k'
	{
		return UVolume2D(x / k, y / k);
	}
	friend std::ostream &operator<<(std::ostream &output, const UVolume2D &temp)//friend��<<���������Ϊ��Ԫ����
	{

		output << "( " << temp.x << " , " << temp.y << " ) ";
		return output;
	}
};

#endif
