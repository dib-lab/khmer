function createSourceLinks() {
    $('.method_details_list .source_code').
        before("<span class='showSource'>[<a href='#' class='toggleSource'>View source</a>]</span>");
    $('.toggleSource').toggle(function() {
       $(this).parent().nextAll('.source_code').slideDown(100);
       $(this).text("Hide source");
    },
    function() {
        $(this).parent().nextAll('.source_code').slideUp(100);
        $(this).text("View source");
    });
}

function createDefineLinks() {
    var tHeight = 0;
    $('.defines').after(" <a href='#' class='toggleDefines'>more...</a>");
    $('.toggleDefines').toggle(function() {
        tHeight = $(this).parent().prev().height();
        $(this).prev().show();
        $(this).parent().prev().height($(this).parent().height());
        $(this).text("(less)");
    },
    function() {
        $(this).prev().hide();
        $(this).parent().prev().height(tHeight);
        $(this).text("more...");
    });
}

function createFullTreeLinks() {
    var tHeight = 0;
    $('.inheritanceTree').toggle(function() {
        tHeight = $(this).parent().prev().height();
        $(this).parent().toggleClass('showAll');
        $(this).text("(hide)");
        $(this).parent().prev().height($(this).parent().height());
    },
    function() {
        $(this).parent().toggleClass('showAll');
        $(this).parent().prev().height(tHeight);
        $(this).text("show all");
    });
}

function searchFrameLinks() {
  $('.full_list_link').click(function() {
    toggleSearchFrame(this, $(this).attr('href'));
    return false;
  });
}

function toggleSearchFrame(id, link) {
  var frame = $('#search_frame');
  $('#search a').removeClass('active').addClass('inactive');
  if (frame.attr('src') == link && frame.css('display') != "none") {
    frame.slideUp(100);
    $('#search a').removeClass('active inactive');
  }
  else {
    $(id).addClass('active').removeClass('inactive');
    frame.attr('src', link).slideDown(100);
  }
}

function framesInit() {
  if (hasFrames) {
    document.body.className += ' frames';
    $('#menu .noframes a').attr('href', document.location);
    try {
        window.top.document.title = $('html head title').text();
    } catch(e) {
        // some browsers like Chrome don't allow this cross-frame access if using file://
    }
  }
  else {
    $('#menu .noframes a').text('frames').attr('href', framesUrl);
  }
}

function fixOutsideWorldLinks() {
  $('a').each(function() {
    if (window.location.host != this.host) this.target = '_parent';
  });
}

function generateTOC() {
  var _toc = $('<ol class="nav"><li class="meta"><a href="#top" class="top"><i class="fa fa-chevron-circle-up"></i> Top</a><a href="page_mainpage.html" class="home"><i class="fa fa-home"></i> Home</a></li></ol>');
  
  var show = false;
  var toc = _toc;
  var counter = 0;
  var tags = ['h2', 'h3', 'h4', 'h5', 'h6'];
  var i;
  if ($('h1').length > 1) tags.unshift('h1');
  for (i = 0; i < tags.length; i++) { tags[i] = tags[i]; }
  var lastTag = parseInt(tags[0][1], 10);

  // iterates through all relevant Hx tags
  $(tags.join(', ')).each(function() {
    if ($(this).parents('.method_details .docstring').length != 0) return;
    if (this.id == "filecontents") return;
    if ($(this).parents('.modal').length != 0) return;
    
    show = true;
    var thisTag = parseInt(this.tagName[1], 10);
    if (this.id.length === 0) {
      var proposedId = $(this).attr('toc-id');
      if (typeof(proposedId) != "undefined") this.id = proposedId;
      else {
        var proposedId = $(this).text().replace(/[^a-z0-9-]/ig, '_');
        if ($('#' + proposedId).length > 0) { proposedId += counter; counter++; }
        this.id = proposedId;
      }
    }
    if (thisTag > lastTag) {
      for (i = 0; i < thisTag - lastTag; i++) {
        var tmp = $('<ol/>'); toc.append(tmp); toc = tmp;
      }
    }
    if (thisTag < lastTag) {
      for (i = 0; i < lastTag - thisTag; i++) toc = toc.parent();
    }
    var title = $(this).attr('toc-title');
    if (typeof(title) == "undefined") title = $(this).text();
    toc.append('<li data-toc="' + ($(this).data('toc') || 'show') + '"><a href="#' + this.id + '">' + title + '</a></li>');
    lastTag = thisTag;
  });
  if (!show) return;
  var html = '<div id="toc"><p class="title"><strong>Table of Contents</strong></p></div>';
  if($('h1').length > 0) {
  	$('h1').first().after(html);
  } else {
  	$('#content').prepend(html);
  }
  
  // hides all items to be hidden (determined by their data-toc attribute) and the optionally following list of sub items
  $(_toc).find('[data-toc=hidden]').each(function() {
  	var $hiddenNavItem = $(this);
  	var $parentList = $hiddenNavItem.parent('ol, ul');
  	if($hiddenNavItem, $hiddenNavItem.next().prop("tagName") == 'OL') $hiddenNavItem.next().remove();
  	$hiddenNavItem.remove();
  	if($parentList.children().length == 0) $parentList.remove();
  });
  $('#toc').append(_toc);
}

$(framesInit);
$(createSourceLinks);
$(createDefineLinks);
$(createFullTreeLinks);
$(searchFrameLinks);
$(fixOutsideWorldLinks);
$(generateTOC);
