/**
 * @file SDL_plot.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements a 2D plot window by using SDL2
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

#include <iostream>

#include "SDL_plot.hpp"
#include "roboto_regular.hpp"
#include "clones_icon.hpp"

namespace CLONES
{

namespace UI
{

SDLWindow::SDLWindow(const unsigned int width, const unsigned int height, const std::string& name):
	Plot2DWindow(width,height,name)
{
	if(SDL_Init(SDL_INIT_VIDEO) < 0) {
		throw std::runtime_error( "SDL initialization error" );
	}

	if (TTF_Init() < 0 ) {
		throw std::runtime_error("SDL_ttf initialization error");
	}

	SDL_RWops* fontMem = SDL_RWFromConstMem(roboto_regular_ttf, sizeof(roboto_regular_ttf));
	if ( !fontMem ) {
		throw std::runtime_error("Failed to load font.");
	}

	font = TTF_OpenFontRW(fontMem, 1, 24);
	if ( !font ) {
		throw std::runtime_error("Failed to load font.");
	}

	SDL_SetHint( SDL_HINT_RENDER_SCALE_QUALITY, "1");

	//Create window
	window = SDL_CreateWindow(name.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
								width, height, SDL_WINDOW_SHOWN);
	if (window == NULL) {
		throw std::runtime_error( "Window could not be created" );
	}

    SDL_RWops *iconMem = SDL_RWFromMem(CLONES_icon, sizeof(CLONES_icon));
	SDL_Surface* icon = IMG_Load_RW(iconMem, 1);

	SDL_SetWindowIcon(window, icon);

	//Create renderer for window
	renderer = SDL_CreateRenderer( window, -1, SDL_RENDERER_ACCELERATED );
	if( renderer == NULL ) {
		throw std::runtime_error( "Renderer could not be created" );
	}
	//Initialize renderer color
	use_color(this->background_color);
}


void SDLWindow::draw_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
					           const unsigned int width, const unsigned int height, const unsigned int thickness)
{
	auto top_left_x = static_cast<int>(upper_left_x)-static_cast<int>(thickness);
	auto top_left_y = static_cast<int>(upper_left_y)-static_cast<int>(thickness);
	auto lower_left_x = top_left_x+static_cast<int>(width+thickness);
	auto lower_left_y = top_left_y+static_cast<int>(height+thickness);

	SDL_Rect outer  = { top_left_x, top_left_y,
					    static_cast<int>(width+2*thickness),
					    static_cast<int>(thickness) };

	SDL_RenderFillRect( renderer, &outer );

	outer  = { top_left_x, top_left_y,
			   static_cast<int>(thickness),
			   static_cast<int>(height+2*thickness) };
			
	SDL_RenderFillRect( renderer, &outer );

	outer  = { top_left_x, lower_left_y,
			   static_cast<int>(width+2*thickness),
			   static_cast<int>(thickness) };

	SDL_RenderFillRect( renderer, &outer );

	outer  = { lower_left_x, top_left_y,
			   static_cast<int>(thickness),
			   static_cast<int>(height+2*thickness) };

	SDL_RenderFillRect( renderer, &outer );
}

void SDLWindow::draw_filled_rectangle(const unsigned int upper_left_x, const unsigned int upper_left_y,
							          const unsigned int width, const unsigned int height)
{
	SDL_Rect rect  = { static_cast<int>(upper_left_x),
						static_cast<int>(upper_left_y),
					    static_cast<int>(width),
					    static_cast<int>(height) };

	SDL_RenderFillRect( renderer, &rect );
}

void SDLWindow::delete_point(const unsigned int x, const unsigned int y)
{
	Color color = get_color();
	
	// set background color
	use_color(this->background_color);

	// plot the point
	draw_point(x,y);

	use_color(color);
}

void SDLWindow::get_text_size(const std::string& text, unsigned int& width, unsigned int& height)
{
	int w, h;
	TTF_SizeText(font, text.c_str(), &w, &h);

	width = static_cast<unsigned int>(w);
	height = static_cast<unsigned int>(h);
}

void SDLWindow::draw_text(const std::string& text,
						  const unsigned int upper_left_x, const unsigned int upper_left_y)
{
	Color color = get_color();

	SDL_Surface* msg_sfc = TTF_RenderText_Blended(font, text.c_str(),
												  {color.red, color.green,
												   color.blue, color.alpha});

	// now you can convert it into a texture
	SDL_Texture* msg = SDL_CreateTextureFromSurface(renderer, msg_sfc);

	unsigned int width, height;

	get_text_size(text, width, height);

	SDL_Rect msg_rect  = { static_cast<int>(upper_left_x),
						   static_cast<int>(upper_left_y),
					       static_cast<int>(width),
					       static_cast<int>(height) };

	SDL_RenderCopy(renderer, msg, NULL, &msg_rect);

	SDL_FreeSurface(msg_sfc);
	SDL_DestroyTexture(msg);
}

void SDLWindow::update()
{
	SDL_RenderPresent( renderer );

	SDL_Event e;
	while(SDL_PollEvent(&e)!=0) {
		//User requests quit
		if( e.type == SDL_QUIT ) {
			this->w_closed = true;

			return;
		}
    }
}

SDLWindow::~SDLWindow()
{
	TTF_CloseFont(font);
	
	SDL_DestroyRenderer( renderer );
	SDL_DestroyWindow( window );

	TTF_Quit();
	IMG_Quit();
	SDL_Quit();
}

}

}