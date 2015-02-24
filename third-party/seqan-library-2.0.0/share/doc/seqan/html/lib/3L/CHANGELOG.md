# CHANGELOG

## 1.4.4
* @font-face has the new parameter - 'filename'. Previously the filename was the same as font-name.
* .hue-rotation() changed to .hue-rotate() since it's the valid name of the filter.
* removed -o-filters since Opera's moving to Webkit.
* the fallback color in gradients is present only if three parameters are used in gradients mixin.
* fallback color in grdients use background-color property instead of background.
* removed Bootstrap integration files.

## 1.4.3-beta
* filter with unlimited arguments.
* specyfic filter functions now done right.
* border-image and its specyfic properties.
* updates in documentation.
* @font-face added.

## 1.4.2-beta (aka Lots of Love for Opera)
* deleted -moz- and -o-perspective.
* added perspective origin.
* deleted -o- in gradients.
* deleted -o- in transforms - Opera either doesn't require prefix for 2D or doesn't support transformations 3D.
* deleted -moz- for 3D transforms.
* deleted -o- transitions.
* multiple specyfic transition properties.
* deleted -o- animations.
* updates and changes in documentation.
* transform-perspective added.

## 1.4.1-beta (aka Lots of Love for Firefox)
* transform style with value preserve-3d is no longer included in 3D transformation and is moved to its own mixin.
* deleted -mox-transition[-*] properties since they're no longer required.
* deleted -mox-animation[-*] properties and -moz-keyframes since they're no longer required.
* deleted -khtml-user-select and -webkit-touch-callout from .user-select() mixin.
* added -moz-border-image for Firefox 3.6.
* deleted -moz- prefixed flex-box properties. Firefox supports the new syntax since version 22 (Aurora channel as for today).

## 1.4.0-beta
* all text blocks in 3L are now preceded by // so they won't appear in a compiled CSS even if you won't minimize your file.
* bugfix for IE gradients (thanks to [Daniel Garcia](https://github.com/zlapper)).
* added animation's and transition's specyfic properties.
* added hasLayout().
* added .selection().
* added .backface-visibility().
* added .appearance() mixin.
* added flex-box's mixins.
* added .user-select() mixin - made by [Daniel Garcia](https://github.com/zlapper).
* added .border-image().
* multi-purpose .gradient() mixin [BETA].
* added filter property.
* unlimited box-shadows and transitions.
* skew is a non standard property, co the .skew() mixin has been changed so now when you use one argument the transformation will be skewX, and if you use two arguments it will be skewX and skewY (transform: skewX(), skewY()).
* added option to include normalize.css, reset.css and HTML5 Boiler Plate's default stylesheet.
* reorgnisation of the folder structure.

## 1.3.3
* -ms-transform back to biznis since IE9 requires that.

## 1.3.2
* added unprefixed transition: transform, since we now know it's safe and future-proof.

## 1.3.1
* fixed bug in gradients. Thanks to [Daniel Garcia](https://github.com/zlapper) for pointing that out!

## 1.3.0
* removed -webkit-background-clip, -webkit-background-size and -webkit-gradient. No longer required.
* removed flex-box mixins. The specification changed, but there's no browsers' support for now. Will be back in the future.
* in gradients, background property changed to background-image.
* gradients' syntax updated to the new specification.

## 1.2.1

* removed -ms- prefixed gradients, transform, transition and animation (together with @keyframes) since [they're no longer required](http://blogs.msdn.com/b/ie/archive/2012/05/31/windows-release-preview-the-sixth-ie10-platform-preview.aspx)

## 1.2.0

* new debbuging rules: a[href="#"], div:empty, span:empty, li:empty, p:empty, td:empty, th:empty, *[title=""], *[class=""] and *[id=""].

## 1.1.1

* mixins updated to be compatible with new variadic arguments - (https://gist.github.com/1933613).
* many @arguments changed to semantic names.

## 1.0.1

* .transition-transform mixin now doesn't require more properties to be declared - you can now use if you want to transition only the transform property.