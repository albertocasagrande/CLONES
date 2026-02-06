/**
 * @file SDL_plot.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines a 2D plot window by using SDL2
 * @version 1.1
 * @date 2026-02-06
 *
 * @copyright Copyright (c) 2023-2026
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef __CLONES_SDL_PLOT__
#define __CLONES_SDL_PLOT__

#include <SDL.h>
#include <SDL_image.h>
#include <SDL_ttf.h>

#include <iostream>

#include "plot_2D.hpp"

namespace CLONES
{

namespace UI
{


/**
 * @brief SDL 2D plotting windows
 */
class SDLWindow : public Plot2DWindow
{
	enum {
		RED = 0,
		GREEN = 1,
		BLUE = 2,
		ALPHA
	};

	SDL_Window* window;         //!< the SDL graphical window
	SDL_Renderer* renderer;     //!< the SDL renderer

	TTF_Font* font;				//!< the text font
public:
	/**
	 * @brief The constructor
	 *
	 * @param width is the width of the plot window
	 * @param height is the height of the plot window
	 * @param name is the name of the plot window
	 */
	SDLWindow(const unsigned int width, const unsigned int height, const std::string& name);

	/**
	 * @brief Clear the window
	 */
	inline void clear()
	{
		use_color(this->background_color);
		SDL_RenderClear( renderer );
	}

	/**
	 * @brief Get the selected color
	 *
	 * @return the selected color
	 */
	inline Color get_color() const
	{
		Color color;
		SDL_GetRenderDrawColor( renderer, &color.red, &color.green,
										&color.blue, &color.alpha );

		return color;
	}

	/**
	 * @brief Use a color
	 *
	 * @param color is the color to be used
	 */
	inline void use_color(const Color& color)
	{
		SDL_SetRenderDrawColor( renderer, color.red, color.green, color.blue, color.alpha );
	}

	/**
	 * @brief Draw a point
	 *
	 * @param x is the x-axis position of the point
	 * @param y is the y-axis position of the point
	 */
	inline void draw_point(const unsigned int x, const unsigned int y)
	{
		SDL_RenderDrawPoint( renderer, x, y );
	}

	/**
	 * @brief Draw a rectangle
	 *
	 * @param upper_left_x is the x-axis position of the rectangle upper left corner
	 * @param upper_left_y is the y-axis position of the rectangle upper left corner
	 * @param width is the rectangle width
	 * @param height is the rectangle height
	 * @param thickness is the rectangle thickness
	 */
	void draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
						const unsigned int width, const unsigned int height, const unsigned int thickness=1);

	/**
	 * @brief Draw a filled rectangle
	 *
	 * @param upper_left_x is the x-axis position of the rectangle upper left corner
	 * @param upper_left_y is the y-axis position of the rectangle upper left corner
	 * @param width is the rectangle width
	 * @param height is the rectangle height
	 */
	void draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
							   const unsigned int width, const unsigned int height);

	/**
	 * @brief Get the size of the a string
	 *
	 * @param[in] text is the string to be drawn
	 * @param[out] width is the width of the string
	 * @param[out] height is the height of the string
	 */
	void get_text_size(const std::string& text, unsigned int& width, unsigned int& height);

	/**
	 * @brief Draw a string
	 *
	 * @param text is the string to be drawn
	 * @param upper_left_x is the x-axis position of the string upper left corner
	 * @param upper_left_y is the y-axis position of the string upper left corner
	 */
	void draw_text(const std::string& text, const unsigned int upper_left_x, const unsigned int upper_left_y);

	/**
	 * @brief Delete a point
	 *
	 * This method delete a point in the plot by drawing over it a point
	 * whose color is the set background color.
	 *
	 * @param x is the x-axis position of the point
	 * @param y is the y-axis position of the point
	 */
	void delete_point(const unsigned int x, const unsigned int y);

	/**
	 * @brief Update the plot in the window
	 */
	void update();

	/**
	 * @brief Test whenever the plotting window is waiting for any event
	 *
	 * @return `true` if and only if the plotting window has not been closed
	 */
	inline bool waiting_end() const
	{
		return !this->closed();
	}

	/**
	 * @brief The SDLWindow destroyer
	 */
	~SDLWindow();
};

}	// UI

}	// CLONES

#endif // __CLONES_SDL_PLOT__