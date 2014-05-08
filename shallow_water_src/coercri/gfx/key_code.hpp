/*
 * FILE:
 *   key_code.hpp
 *
 * PURPOSE:
 *   Key events, keycodes, modifier codes, function to convert keycode
 *   to key name.
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

#ifndef COERCRI_KEY_CODE_HPP
#define COERCRI_KEY_CODE_HPP

#include <string>

namespace Coercri {

    // Raw key codes
    
    // NOTE: New RawKeys should be added at the END, so that programs
    // relying on the existing numeric codes do not break.
    
    enum RawKey {
        RK_0,
        RK_1,
        RK_2,
        RK_3,
        RK_4,
        RK_5,
        RK_6,
        RK_7,
        RK_8,
        RK_9,
        RK_A,
        RK_AMPERSAND,
        RK_ASTERISK,
        RK_AT,
        RK_B,
        RK_BACKQUOTE,
        RK_BACKSLASH,
        RK_BACKSPACE,
        RK_BREAK,
        RK_C,
        RK_CAPS_LOCK,
        RK_CARET,
        RK_CLEAR,
        RK_COLON,
        RK_COMMA,
        RK_COMPOSE,
        RK_D,
        RK_DELETE,
        RK_DOLLAR,
        RK_DOUBLE_QUOTE,
        RK_DOWN,
        RK_E,
        RK_END,
        RK_EQUALS,
        RK_ESCAPE,
        RK_EURO,
        RK_EXCLAIM,
        RK_F,
        RK_F1,
        RK_F2,
        RK_F3,
        RK_F4,
        RK_F5,
        RK_F6,
        RK_F7,
        RK_F8,
        RK_F9,
        RK_F10,
        RK_F11,
        RK_F12,
        RK_F13,
        RK_F14,
        RK_F15,
        RK_G,
        RK_GREATER,
        RK_H,
        RK_HASH,
        RK_HELP,
        RK_HOME,
        RK_I,
        RK_INSERT,
        RK_J,
        RK_K,
        RK_KP_0,   // first numeric keypad key.
        RK_KP_1,
        RK_KP_2,
        RK_KP_3,
        RK_KP_4,
        RK_KP_5,
        RK_KP_6,
        RK_KP_7,
        RK_KP_8,
        RK_KP_9,
        RK_KP_DIVIDE,
        RK_KP_ENTER,
        RK_KP_EQUALS,
        RK_KP_MINUS,
        RK_KP_MULTIPLY,
        RK_KP_PERIOD,
        RK_KP_PLUS,   // last numeric keypad key.
        RK_L,
        RK_LEFT,
        RK_LEFT_ALT,
        RK_LEFT_BRACKET,
        RK_LEFT_CONTROL,
        RK_LEFT_META,
        RK_LEFT_PAREN,
        RK_LEFT_SHIFT,
        RK_LEFT_WINDOWS,
        RK_LESS,
        RK_M,
        RK_MENU,
        RK_MODE,
        RK_MINUS,
        RK_N,
        RK_NUM_LOCK,
        RK_O,
        RK_P,
        RK_PAGE_DOWN,
        RK_PAGE_UP,
        RK_PAUSE,
        RK_PERIOD,
        RK_PLUS,
        RK_POWER,
        RK_PRINT_SCREEN,
        RK_Q,
        RK_QUESTION,
        RK_R,
        RK_RETURN,
        RK_RIGHT,
        RK_RIGHT_ALT,
        RK_RIGHT_BRACKET,
        RK_RIGHT_CONTROL,
        RK_RIGHT_META,
        RK_RIGHT_PAREN,
        RK_RIGHT_SHIFT,
        RK_RIGHT_WINDOWS,
        RK_S,
        RK_SCROLL_LOCK,
        RK_SEMICOLON,
        RK_SINGLE_QUOTE,
        RK_SLASH,
        RK_SPACE,
        RK_SYSREQ,
        RK_T,
        RK_TAB,
        RK_U,
        RK_UNDERSCORE,
        RK_UNDO,
        RK_UP,
        RK_V,
        RK_W,
        RK_X,
        RK_Y,
        RK_Z,
        RK_UNKNOWN
    };

    std::string RawKeyName(RawKey rk);


    // Cooked key codes
    enum CookedKey {
        CK_CHARACTER,

        CK_BACKSPACE,
        CK_DELETE,
        CK_DOWN,
        CK_END,
        CK_ESCAPE,
        CK_F1,
        CK_F2,
        CK_F3,
        CK_F4,
        CK_F5,
        CK_F6,
        CK_F7,
        CK_F8,
        CK_F9,
        CK_F10,
        CK_F11,
        CK_F12,
        CK_F13,
        CK_F14,
        CK_F15,
        CK_HELP,
        CK_HOME,
        CK_INSERT,
        CK_LEFT,
        CK_LEFT_WINDOWS,
        CK_MENU,
        CK_MODE,
        CK_PAGE_DOWN,
        CK_PAGE_UP,
        CK_PAUSE,
        CK_PRINT_SCREEN,
        CK_RETURN,
        CK_RIGHT,
        CK_RIGHT_WINDOWS,
        CK_TAB,
        CK_UP,

        CK_UNKNOWN
    };

    // Modifiers for cooked keys
    // (OR'd together)
    enum KeyModifier {
        KM_ALT = 1,
        KM_CONTROL = 2,
        KM_SHIFT = 4
    };
}

#endif
