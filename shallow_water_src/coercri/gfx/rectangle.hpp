/*
 * FILE:
 *   rectangle.hpp
 *
 * PURPOSE:
 *   Rectangle class
 *
 * AUTHOR:
 *   Stephen Thompson <stephen@solarflare.org.uk>
 *
 * COPYRIGHT:
 *   Copyright (C) Stephen Thompson, 2008 - 2009.
 *
 *   This file is part of the "Coercri" software library. Usage of "Coercri"
 *   is permitted under the terms of the Boost Software License, Version 1.0, 
 *   the text of which is displayed below.
 *
 *   Boost Software License - Version 1.0 - August 17th, 2003
 *
 *   Permission is hereby granted, free of charge, to any person or organization
 *   obtaining a copy of the software and accompanying documentation covered by
 *   this license (the "Software") to use, reproduce, display, distribute,
 *   execute, and transmit the Software, and to prepare derivative works of the
 *   Software, and to permit third-parties to whom the Software is furnished to
 *   do so, all subject to the following:
 *
 *   The copyright notices in the Software and this entire statement, including
 *   the above license grant, this restriction and the following disclaimer,
 *   must be included in all copies of the Software, in whole or in part, and
 *   all derivative works of the Software, unless such copies or derivative
 *   works are solely in the form of machine-executable object code generated by
 *   a source language processor.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *   FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 *   SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 *   FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 *   DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef COERCRI_RECTANGLE_HPP
#define COERCRI_RECTANGLE_HPP

#include <algorithm>

namespace Coercri {

    class Rectangle {
    public:
        // Note invariant: width >= 0 && height >= 0

        // Constructors
        Rectangle()
            : left(0), top(0), width(0), height(0) { }
        Rectangle(int l, int t, int w, int h)
        : left(l), top(t), width(w), height(h) {
            if (width < 0) width = 0; if (height < 0) height = 0;
        }

        // query positions
        int getLeft() const { return left; }
        int getRight() const { return left + width; }
        int getTop() const { return top; }
        int getBottom() const { return top + height; }

        // query size
        int getWidth() const { return width; }
        int getHeight() const { return height; }

        // "degenerate" means zero area.
        bool isDegenerate() const { return width == 0 || height == 0; }

        // Move/Resize
        void moveTo(int l, int t) {
            left = l; top = t;
        }
        void translate(int x, int y) {
            left += x; top += y;
        }
        void setSize(int w, int h) {  // fixes top left corner
            width = w; if (width < 0) width = 0;
            height = h; if (height < 0) height = 0;
        }
        void addToSize(int w, int h) {
            width += w; if (width < 0) width = 0;
            height += h; if (height < 0) height = 0;
        }
        
    private:
        int left, top;
        int width, height;
    };

    // Intersect two rectangles
    inline Rectangle IntersectRects(const Rectangle &r1, const Rectangle &r2)
    {
        const int new_left = std::max(r1.getLeft(), r2.getLeft());
        const int new_right = std::min(r1.getRight(), r2.getRight());
        const int new_top = std::max(r1.getTop(), r2.getTop());
        const int new_bottom = std::min(r1.getBottom(), r2.getBottom());
        return Rectangle(new_left, new_top, new_right - new_left, new_bottom - new_top);
    }

    // Check if a point is in a rectangle
    inline bool PointInRect(const Rectangle &r, int x, int y)
    {
        return x >= r.getLeft() && x < r.getRight() && y >= r.getTop() && y < r.getBottom();
    }
}

#endif
