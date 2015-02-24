# 3L was made for YOU to help you create awesome websites and fill the Internet with excessive amount of Love! ♥

* Author: [Mateusz Kocz](http://radiatingstar.com)
* [Watch 3L on Github](https://github.com/mateuszkocz/3l)
* [Submit a bug issue](https://github.com/mateuszkocz/3l/issues?state=open)
* [Download .zip](https://github.com/mateuszkocz/3l/archive/master.zip)

3L is MIT licensed.
Includes normalize.css, reset.css and some code from HTML5 Boilerplate.
For the licencess refer to the LICENCES.md.

3L version: 1.4.4 (2013.06.03)
LESS version: 1.3.3
Reset.css version: 2.0
Normalize.css version: 2.1.2
HTML5BP's CSS version: 4.0.1

[Get your own LESS.js](http://lesscss.org/)

---

# HOW TO USE

## Basic use
[Download 3L.zip file](https://github.com/mateuszkocz/3l/archive/master.zip), unzip it and place it in your project. You can start editing `style.less` file or `@import 3L.less` into your previously created LESS stylesheet. Use `.normalize()` or `.reset()` classes if you want.

## Namespacing
If you're using anoter mixins library you might want to put 3L into a namespace so the two libraries won't clash. Just type `#3L {@import '3L/3L';}` and the 3L will be fully contained in its own namespace. Access the mixins with `#3L > .mixin()`.


## Animations
`@import "animation[1-5]"` (they're in 3L/assets/animationsins) to your stylesheet, create a class `.animation[1-5]() {/* @keyframes properties */}` and declare your animation. 3L will make the prefixed @keyframes for you. Now just use this animation in any element you want with `.animation()` mixin.

## Refer to the documentation
The 3L.less file has all documentation included. If you want to know more how to use a mixin, what parameters it takes, what browsers it supports or where you could find more information about a CSS property just find the corresponding section in the file or refer to the [3L wiki on GitHub](https://github.com/mateuszkocz/3l/wiki).

## Compile
All .less files should be compiled to .css. You can use native LESS compiler or you can try [Winless](http://winless.org/), [Prepross](http://alphapixels.com/prepros/) and [CodeKit](http://incident57.com/codekit/). Your output CSS file will be clean. Only the used mixins from 3L will be included, nothing else.

---

# ANIMATIONS

Using animations with 3L is pretty easy. Copy the animations less files (animationX.less) to the folder where you have your style sheet. In your code just type `@import "animation1"` (or any other number), create a class `.animation1` (or any other) and declare your animation. 3L will make @keyframes for you. Now just use this animation in any element you want with `.animation()` class.
```
	@import "animation1";
	.animation1 () {
		/* your @keyframes rules */
		}
	.someClassName {
		.animation(.animation1 1s);
		}
```
---

# FAQ and TROUBLESHOOTING

Compiler I use fail to compile 3L!

Unfortunately not all compilers can deal with some of the 3L classes. The hardest to compile are opacity in percentages (it uses a small piece of JavaSript), keyframes and guarded mixins. I strongly encourage you to use WinLess (also the online version) as it works very well with all 3L classes.

Other solution would be to delete the option to declare opacity in percentages (it's the biggest class in opacity block in 3L.less file and uses JavaScript). If it still doesn't help, you seriously need to consider getting a better compiler compatible with all LESS functionality!
