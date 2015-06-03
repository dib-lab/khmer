/**
 * Various
 */
(function ($) {

	$(document).ready(function () {

	    // hightlight nav items based on scroll area
	    $('body').scrollspy({ target: '#toc', offset: 50 });
	    
	    var id;
	    $(window).resize(function() {
            clearTimeout(id);
            id = setTimeout(function() { $('body').scrollspy('refresh'); }, 500);
        });
        
        // shows 'open in frameset' link if opened separately
		if(window == window.parent && window.name != 'list') {
        	$('#content').prepend('<div class="open-in-frame alert alert-info"><a href="index.html?p=' + $('html').attr('data-page') + window.location.hash + '"><strong>Looking for a different entry?</strong> Unhide the navigation bar and start your search.</a></div>'); 
        }
        
        // if loaded in a frame, checks the URI's p parameter and uses it to load the
        // specified page in the right frame
        // (e.g. docs.seqan.de/index.html?p=String#Example will open the example section of the class String in the main frame) 
        if(window != window.parent && window.name == 'list') {
            try {
                var redirectTo = null;
                var hash = $.urlHash(window.parent.location);
                if($.urlParam('p', window.parent.location)) {
                    var p = $.urlParam('p', window.parent.location).split('/')[0];
                    if (p.indexOf('::') != -1)
                    {
                        var tmp = p;
                        p = tmp.split('::')[0];
                        hash = '::' + tmp.split('::')[1];
                        console.log('p == ' + p + ' -- hash = ' + hash);
                    }
                    if(window.lookup.hasOwnProperty(p)) {
                        redirectTo = window.lookup[p] + '.html#' + encodeURIComponent(p + hash);
                    } else {
                        $(window.parent['main'].document).find('#content').prepend('<div class="open-in-frame alert alert-danger">Could not find page for <strong>' + p + '</strong></div>'); 
                        // TODO: start search using search form for p
                    }
                } else {
                    if(hash.length > 1) {
                        redirectTo = hash.substr(1) + '.html';
                    }
                }

                if(redirectTo) {
                    window.parent['main'].location = redirectTo;
                }
    		} catch(e) {
    		    // some browsers like Chrome don't allow this cross-frame access if using file://
    		}
        }
        
        // adds a close link to the search/list frame
		if(window != window.parent && window.name == 'list') {
			try {
        		$('<button class="close" style="position: absolute; top: 5px; right: 5px; line-height: .7em;">&times;</button>').prependTo('body').click(function() {
        			window.parent.location = window.parent['main'].location;
        		});
    		} catch(e) {
    		    // some browsers like Chrome don't allow this cross-frame access if using file://
    		}
        }
 
        // tooltips
        $('[title]:not([href])').tooltip({ container: 'body' });

        // smooth scrolling
        //$('a[href*=#]:not([href=#])').smoothScroll({ offset: -20 });
        
        // autofocus search field
        if($('html').hasClass('list')) {
        	window.setTimeout(function() {
        		$('input[type=search]').focus();
        	}, 50);
		}
		
    });

})(jQuery);

/**
 * Get URL parameter functionality
 */
(function ($) {
	$.extend({
		urlParam: function(name, location) {
			if(!location) location = window.location;
            return decodeURIComponent((new RegExp('[?|&]' + name + '=' + '([^&;]+?)(&|#|;|$)').exec(location.search) || [, ""])[1].replace(/\+/g, '%20')) || null;
		},
		urlHash: function(location) {
			if(!location) location = window.location + '';
			else location += '';
			
			var index = location.indexOf('#');
            if(index >= 0) {
            	return location.substr(index);
            } else {
            	return '';
            }
		}
	});
})(jQuery);

/**
 * Language Entity Labels
 */
(function ($) {

	$.fn.extend({
		createLangEntityLabel: function(langEntity) {
            if(!window.langEntities) return;
            var entry = window.langEntities[langEntity];
			if (!entry) entry = window.langEntities['unknown'];
			
			if(langEntity == 'tutorial') return $('<span>' + entry.ideogram + '</span>');
			return $('<a href="page_LanguageEntities.html#' + langEntity + '">' + entry.ideogram + '</a>');
		},
		
		pimpLangEntityLabel: function(langEntity) {
    		$this = $(this);
    		if(jQuery.inArray($this.prop('tagName'), ['A']) != -1) {
    		    // a tags may and should not be nested
			    $this.wrap('<span data-lang-entity="' + langEntity + '" data-pimped="true"/>')
					.removeAttr('data-lang-entity')
					.before($().createLangEntityLabel(langEntity));
			} else {
    			$this.wrapInner('<span/>').prepend($().createLangEntityLabel(langEntity));
			}
		},
		
		/**
		 * Annotates all tags with a data-lang-entity(-container) attribute the following way:
		 * 1) if the tag has the data-lang-entity attribute,
		 *      it will be prefixed with a colored lang entity label
		 * 2) if the tag has the data-lang-entity-container attribute,
         *      the highest headline level will be prefixed,
         *      whereas all other headlines won't be prefixed (by removing their eventually set data-lang-entity attribute)
		 */
		pimpLangEntityLabels: function() {
			return this.each(function() {
				$(this).find('[data-lang-entity]').each(function () {
					var $this = $(this);
					if($this.attr('data-pimped')) return true;
					//if($this.parents('[data-lang-entity-container]').length > 0) return true;
            
					var langEntity = $this.attr('data-lang-entity');
					$this.pimpLangEntityLabel(langEntity);				
				});
			}).each(function() {
			/*
			    $(this).find('[data-lang-entity-container]').each(function () {
    			    // only use one headline level
    			    var headlineHandled = false;
    			    for(var i=1; i<=6; i++) {
        			    $el = $(this).find('h' + i);
        			    if(!headlineHandled) {
            			    if($el.length > 0) {
            					if($el.attr('data-pimped')) return true;
            					
            					var langEntity = $el.parents('[data-lang-entity-container]').attr('data-lang-entity-container');
            					$el.attr('data-lang-entity', langEntity);
            					$el.pimpLangEntityLabel(langEntity);
            					
                			    headlineHandled = true;
            			    }
        			    } else {
            			    //$el.removeAttr('data-lang-entity');
        			    }
    			    }
    			 });
    			 */
			});
		}
	});

    $(document).ready(function () {
        $('body').on('mouseover', '[data-lang-entity] > :first-child', function() {
            var langEntity = $(this).attr('data-lang-entity') || $(this).parent().attr('data-lang-entity');
            showPopOver(this, langEntity, true);
        });
        
        $('#search').on('mouseover', '[data-lang-entity-container] label', function() {
            var langEntity = $(this).parents('[data-lang-entity-container]').attr('data-lang-entity-container');
            showPopOver(this, langEntity, false);
        });
        
        function showPopOver(el, langEntity, showMore) {
    		var $this = $(el);
    		if($this.attr('data-original-title')) return; // already set up

    		var langEntityData = window.langEntities[langEntity];
    		
    		$this.popover({
				html: true,
				trigger: 'hover',
				template: '<div class="popover" data-lang-entity-container="' + langEntity + '"><div class="arrow"></div><h3 class="popover-title"></h3><div class="popover-content"></div></div>',
				title: langEntityData.name,
				content: function() {
					return '<div class="description">' + langEntityData.description + '</div>' + (showMore && langEntity != 'tutorial' ? '<p class="more">Click now for more information...</p>' : '');
				},
				container: 'body',
				placement: function(tip, element) {
                    var $thisement, above, actualHeight, actualWidth, below, boundBottom, boundLeft, boundRight, boundTop, elementAbove, elementBelow, elementLeft, elementRight, isWithinBounds, left, pos, right;
                    isWithinBounds = function(elementPosition) {
                    return boundTop < elementPosition.top && boundLeft < elementPosition.left && boundRight > (elementPosition.left + actualWidth) && boundBottom > (elementPosition.top + actualHeight);
                    };
                    $element = $(element);
                    pos = $.extend({}, $element.offset(), {
                        width: element.offsetWidth,
                        height: element.offsetHeight
                    });
                    actualWidth = 283;
                    actualHeight = 117;
                    boundTop = $(document).scrollTop() + 40; // DIRTY: takes the small data-lang-entity window "hat" into account
                    boundLeft = $(document).scrollLeft();
                    boundRight = boundLeft + $(window).width();
                    boundBottom = boundTop + $(window).height();
                    
                    elementAbove = {
                        top: pos.top - actualHeight,
                        left: pos.left + pos.width / 2 - actualWidth / 2
                    };
                    elementBelow = {
                        top: pos.top + pos.height,
                        left: pos.left + pos.width / 2 - actualWidth / 2
                    };
                    elementLeft = {
                        top: pos.top + pos.height / 2 - actualHeight / 2,
                        left: pos.left - actualWidth
                    };
                    elementRight = {
                        top: pos.top + pos.height / 2 - actualHeight / 2,
                        left: pos.left + pos.width
                    };
                    above = isWithinBounds(elementAbove);
                    below = isWithinBounds(elementBelow);
                    left = isWithinBounds(elementLeft);
                    right = isWithinBounds(elementRight);
                    if (above) {
                        return "top";
                    } else {
                        if (below) {
                            return "bottom";
                        } else {
                            if (left) {
                                return "left";
                            } else {
                                if (right) {
                                    return "right";
                                } else {
                                    return "bottom";
                                }
                            }
                        }
                    }
                }
			}).popover('show');
    	}
    	
    	if(!$('html').hasClass('list')) {
    		$('html').pimpLangEntityLabels();
    	}
    });

})(jQuery);





/**
 * Permalink Modal
 */
(function ($) {
    $.fn.extend({
		permalinkModal: function(options) {
		
			var settings = $.extend({
				modalId: 'permalinkModal',
				modalTemplate: '<div class="modal fade" id="permalinkModal" tabindex="-1" role="dialog" aria-hidden="true">\
									<div class="modal-dialog">\
										<div class="modal-content">\
											<div class="modal-header">\
												<button type="button" class="close" data-dismiss="modal" aria-hidden="true">&times;</button>\
												<h4 class="modal-title">Permalink</h4>\
											</div>\
											<div class="modal-body">\
												<p>The permalinks for <strong></strong> are:</p>\
												<table>\
													<tr><th>Frame-based</th><td></td><td><button type="button" class="btn btn-xs btn-primary" data-link="true">Copy</button></td></tr>\
													<tr><th>Frame-less</th><td></td><td><button type="button" class="btn btn-xs btn-primary" data-link="true">Copy</button></td></tr>\
													<tr><th>Dox</th><td></td><td><button type="button" class="btn btn-xs btn-primary" data-link="true">Copy</button></td></tr>\
												</table>\
											</div>\
											<div class="modal-footer">\
												<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>\
											</div>\
										</div>\
									</div>\
								</div>',
				elementSelector: '.modal-body strong',
				linkSelectors: ['.modal-body tr:nth-child(1) td:nth-child(2)', '.modal-body tr:nth-child(2) td:nth-child(2)', '.modal-body tr:nth-child(3) td:nth-child(2)'],
				element: 'myElement',
				links: ['#', '#', '#']
			}, options);
					
			function createModal(parent) {
				if($('#' + settings.modalId).length == 0) {
					var $modal = $(settings.modalTemplate)
						.attr('id', settings.modalId)
						.attr('aria-labelledby', settings.modalId + 'Label');
					$modal
						.find('h4')
						.first()
						.attr('id', settings.modalId + 'Label');
					$modal.appendTo(parent);
					
					$modal.find('[data-link]').each(function() {
						var clip = new ZeroClipboard(this, {
							moviePath: "lib/ZeroClipboard/ZeroClipboard.swf"
						});
					
						clip.on("load", function(client) {
							client.on( "complete", function(client, args) {
								$modal.modal('hide');
							});
						});
					});
				}
			}
			
			function changeElement(element) {
				$('#' + settings.modalId + ' ' + settings.elementSelector).html(element);
			}
			
			function changeLinks(links) {
				for(var i=0; i<settings.linkSelectors.length; i++) {
					var html = (links[i].substring(0,1) != '@')
						? '<a href="' + links[i] + '" target="_blank">' + links[i] + '</a>'
						: links[i];
					$('#' + settings.modalId + ' ' + settings.linkSelectors[i]).html(html);
					$($('#' + settings.modalId + ' [data-link]')[i]).attr('data-clipboard-text', links[i]);
				}			
			}
		
			return this.each(function() {
				createModal(this);
				changeElement(settings.element);
				changeLinks(settings.links);
				$('#' + settings.modalId).modal({});
			});
			
		}
	});

    $(document).ready(function () {
    	if($('html').hasClass('list')) return;
    	
    	function handleClick() {
    		var href = window.location.href;
			var fragmentName = $(this).attr('id') || $(this).attr('name') || null;
			var fragment = fragmentName ? '#' + encodeURIComponent(fragmentName) : '';
			$('body').permalinkModal({
				element: fragmentName ? 'fragment ' + fragmentName : 'this page',
				links: [
					href.substring(0, href.lastIndexOf('/')+1) + '?p=' + $('html').data('page') + fragment,
					href.split('#')[0] + fragment,
					'@link ' + $('html').data('page') + fragment + ' @endlink'
				]
			});
    	}
    	
		function permalinks(activate) {
			if(activate) {				
				$('h1,[id],[name]')
					.filter(function() {
						if($.inArray($(this).attr('id'), ['content', 'toc', 'filecontents', 'devModeWindow']) != -1) return false;
						if($(this).hasClass('global-zeroclipboard-container')) return false;
						if($(this).hasClass('modal')) return false;
						if($(this).parents('.modal').length > 0) return false;
						return true; })
					.addClass('permalink')
					.bind('click', handleClick);
				permalinksVisible = true;
			} else {
				$('.permalink').unbind('click', handleClick);
				$('.permalink').removeClass('permalink');
				permalinksVisible = false;
			}
		}
		
		permalinks($.devMode());
		$(document).bind('devMode', function(e) {
			permalinks(e.active);
		});
    });

})(jQuery);




/**
 * Developer Mode
 * Activated / deactivated by pressing Ctrl + Shift at the same time
 */
(function ($) {
	if($('html').hasClass('list')) return;
	
    $.extend({
		devMode: function() {
			var args = Array.prototype.slice.call(arguments);
			if(args.length == 1) {
				var active = args[0] ? true : false;
				localStorage.setItem('devMode', active ? 'true' : 'false');
				$.event.trigger({
					type: 'devMode',
					active: active,
					time: new Date()
				});
				console.log('developer mode: ' + (active ? 'on' : 'off'));
				if(active && $('#devModeWindow').length == 0) {
					$('<div id="devModeWindow"><strong>Developer Mode is active</strong>\
					   <br>Press <code>Ctrl + Shift</code> to deactivate<br></div>')
					   	.append($('<a href="#">Show dox sources</a>').click(function() { $('#doxSources').modal({}).find('.modal-dialog').css('width', '90%'); }))
						.appendTo('body');
				} else {
					$('#devModeWindow').remove();
				}
			} else {
				return localStorage.getItem('devMode') == 'true' ? true : false;
			}
		}
	});
	
	$(document).ready(function () {
		// trigger devMode event at load
		$.devMode($.devMode());
	});
	
	var lastKeys = 0; // last time ctrl + shift was fired - used to detect double fired events (experienced on linux)
	var lastKeysWindow = 500; // time frame within no further key combination is considered
	$(document).bind('keyup keydown', function(e) {
		if(e.ctrlKey && e.shiftKey && lastKeys + lastKeysWindow < new Date().getTime()) {
			lastKeys = new Date().getTime();
			
			if($.devMode()) $.devMode(false);
			else $.devMode(true);
		}
	});
})(jQuery);
		





/**
 * Code Collapse
 */
(function ($) {
    $.fn.extend({
		codeCollapse: function(options) {
		
			var settings = $.extend({
				maxHeight: 200,
				moreLink: '<a class="more">More ...</a>',
				lessLink: '<a class="less">Less ...</a>',
				tolerance: 50 // number of pixels a container's height may exceed before it becomes collapsed
			}, options);
			
			function createMoreLink(box) {
				var $box = $(box);
				
				$box.height('auto');
				var expandedHeight = $box.outerHeight();
				if(expandedHeight <= settings.maxHeight + settings.tolerance) return;
				
				//var srcPath = $box.parents('[data-src-path]').data('src-path');
				//console.log(srcPath, $('[data-src-path="' + srcPath + '.stdout"]').length);
				
				$box.height(settings.maxHeight).css({ overflow: 'hidden' });
				return $(settings.moreLink).click(function() {
					var $link = $(this);
					$box.animate({'height': expandedHeight }, 400);
					$link.fadeOut(400, function() { $link.replaceWith(createLessLink(box)); });
				});
			}
			
			function createLessLink(box) {
				var $box = $(box);
				
				return $(settings.lessLink).click(function() {
					var $link = $(this);
					$box.animate({'height': settings.maxHeight }, 400);
					$link.fadeOut(400, function() { $link.replaceWith(createMoreLink(box)); });
				});
			}
		
			return this.each(function() {
				$(createMoreLink(this)).insertAfter(this);
			});
		}
	});

    $(document).ready(function () {
        $('[data-src-path] pre, pre[data-src-path]').codeCollapse({
        	maxHeight: 77,
        	moreLink: '<a class="more">...</a>',
        	lessLink: '<a class="less">&nbsp;</a>'
        });
    });

})(jQuery);





/**
 * Search Bar
 */
(function ($) {

    function createFilterableSearch($el) {
        // make filter box look fancy
        $el.find('select').multiselect({
            buttonClass: 'btn btn-primary',
            includeSelectAllOption: true,
            selectAllText: "(Un)check all",
            selectAllValue: 'all',
            dropRight: false,
            buttonText: function (checkedOptions) {
                var options = $(arguments[1][0]).find('option').map(function () {
                    return $(this).val();
                });
                $.each(options, function(i){
                    if(options[i] === 'all') {
                        options.splice(i,1);
                        return false;
                    }
                });
                $.each(checkedOptions, function(i){
                    if(checkedOptions[i] === 'all') {
                        checkedOptions.splice(i,1);
                        return false;
                    }
                });

				var $btn = $el.find('button.multiselect');
				if(checkedOptions.length == options.length) {
					$btn.removeClass('btn-warning');
					$btn.addClass('btn-primary');
				} else {
					$btn.addClass('btn-warning');
					$btn.addClass('btn-primary');
				}

                if (options.length == checkedOptions.length) return '<i class="fa fa-filter"></i><small>all visible</small>';
                else if (options.length - checkedOptions.length == 1) return '<i class="fa fa-filter"></i><small>1 excluded</small>';
                else if (checkedOptions.length == 0) return '<i class="fa fa-filter"></i><small>all excluded</small>';
                else return '<i class="fa fa-filter"></i><small>' + (options.length - checkedOptions.length) + ' excluded</small>';
            },
            onChange: function ($element, checked) {
                //$allLabel = $('.multiselect-container').find('li:first-child label');
                //$allLabel.contents().last()[0].textContent = $allLabel.text()[0] == 'U' ? 'Check all' : 'Uncheck all';
            }
        });
        $el.find('.multiselect-container [value=all]').prop('checked', true);
        
        $el.find('button.multiselect').attr('title', '');

        // copies the options value to <li>'s data-lang-entity-container attribute classes of the parent li element (for easier styling)
        $el.find('.multiselect-container input[value]').each(function () {
            $this = $(this);
            if($this.val() == 'all') return;
            $this.parents('li').attr('data-lang-entity-container', $this.val());
            $this.parents('a').click(function (e) {
                if (e.target == this) {
                    // link and not the label or the input was clicked
                    $(this).find('input').click();
                }
            });
        });

        if($el.jsonsearch) $el.jsonsearch({
            numElementsPerPage: -1,
            target: 'main',
            raw: true,
            showUrl: false,
            minimumLength: 1,
            descriptiveWords: 25,
            highlightTerms: true,
            highlightEveryTerm: true,
            output: $("#results"),
            data: window.searchData,
            stopWords: [], // filtered out of query
            replaceWords: [ // words replaced in the query
                []
            ],
            stemWords: [ // silently adds the stem if the corresponding word was found in the query
                ["javascript", "script"]
            ],
            langEntityGroups: [ // TODO create from window[langEntities] (try console.log(window[langEntities])
                ["grouped_typedef", "typedef"],
                ["global_typedef", "typedef"],
                ["member_typedef", "typedef"],
                ["grouped_tag", "tag"],
                ["global_variable", "variable"],
                ["local_variable", "variable"],
                ["member_variable", "variable"]
            ],
            callback: function($form, $results) {
            	if($form.find('input[type=search]').val().length == 0) {
            		$("html").removeClass('shows-results');
            	    $("#results").fadeOut();
            	    $el.find('.pre-action').slideDown();
                } else {
                	$("html").addClass('shows-results');
                    $("#results").fadeIn();
                    $el.find('.pre-action').slideUp();
                }
            }
        });
    }

    $(document).ready(function () {
        createFilterableSearch($('#search'));
        $('#search').fadeIn();
        
        try {
            // search immediately if query was passed within the url
    		var q = decodeURI((RegExp('q=' + '(.+?)(&|$)').exec(parent.location.search)||[,null])[1]);
    		if(q && q != 'null') {
    			$('#search [type=search]').val(q).change().focus();
    		}
        } catch(e) {
            // some browsers like Chrome don't allow this cross-frame access if using file://
        }
    });
    
    // hide form and results until they are pimped
	$('head').append('<style>#search, #results { display: none; }</style>');

})(jQuery);
