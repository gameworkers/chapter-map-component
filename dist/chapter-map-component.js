(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory(require("react"), require("prop-types"), require("@gameworkers/react-simple-maps"), require("redux-tooltip"), require("redux"), require("next-redux-wrapper"));
	else if(typeof define === 'function' && define.amd)
		define(["react", "prop-types", "@gameworkers/react-simple-maps", "redux-tooltip", "redux", "next-redux-wrapper"], factory);
	else if(typeof exports === 'object')
		exports["ChapterMapComponent"] = factory(require("react"), require("prop-types"), require("@gameworkers/react-simple-maps"), require("redux-tooltip"), require("redux"), require("next-redux-wrapper"));
	else
		root["ChapterMapComponent"] = factory(root["React"], root["PropTypes"], root["ReactSimpleMaps"], root["ReduxTooltip"], root["Redux"], root["NextReduxWrapper"]);
})((typeof self !== "undefined" ? self : this), function(__WEBPACK_EXTERNAL_MODULE__0__, __WEBPACK_EXTERNAL_MODULE__1__, __WEBPACK_EXTERNAL_MODULE__2__, __WEBPACK_EXTERNAL_MODULE__3__, __WEBPACK_EXTERNAL_MODULE__4__, __WEBPACK_EXTERNAL_MODULE__5__) {
return /******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, { enumerable: true, get: getter });
/******/ 		}
/******/ 	};
/******/
/******/ 	// define __esModule on exports
/******/ 	__webpack_require__.r = function(exports) {
/******/ 		if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 			Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 		}
/******/ 		Object.defineProperty(exports, '__esModule', { value: true });
/******/ 	};
/******/
/******/ 	// create a fake namespace object
/******/ 	// mode & 1: value is a module id, require it
/******/ 	// mode & 2: merge all properties of value into the ns
/******/ 	// mode & 4: return value when already ns object
/******/ 	// mode & 8|1: behave like require
/******/ 	__webpack_require__.t = function(value, mode) {
/******/ 		if(mode & 1) value = __webpack_require__(value);
/******/ 		if(mode & 8) return value;
/******/ 		if((mode & 4) && typeof value === 'object' && value && value.__esModule) return value;
/******/ 		var ns = Object.create(null);
/******/ 		__webpack_require__.r(ns);
/******/ 		Object.defineProperty(ns, 'default', { enumerable: true, value: value });
/******/ 		if(mode & 2 && typeof value != 'string') for(var key in value) __webpack_require__.d(ns, key, function(key) { return value[key]; }.bind(null, key));
/******/ 		return ns;
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "/dist";
/******/
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = 6);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__0__;

/***/ }),
/* 1 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__1__;

/***/ }),
/* 2 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__2__;

/***/ }),
/* 3 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__3__;

/***/ }),
/* 4 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__4__;

/***/ }),
/* 5 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__5__;

/***/ }),
/* 6 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
__webpack_require__.r(__webpack_exports__);

// EXTERNAL MODULE: external {"root":"React","commonjs":"react","commonjs2":"react","amd":"react"}
var external_root_React_commonjs_react_commonjs2_react_amd_react_ = __webpack_require__(0);
var external_root_React_commonjs_react_commonjs2_react_amd_react_default = /*#__PURE__*/__webpack_require__.n(external_root_React_commonjs_react_commonjs2_react_amd_react_);

// EXTERNAL MODULE: external {"root":"PropTypes","commonjs":"prop-types","commonjs2":"prop-types","amd":"prop-types"}
var external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_ = __webpack_require__(1);
var external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default = /*#__PURE__*/__webpack_require__.n(external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_);

// EXTERNAL MODULE: external {"root":"ReactSimpleMaps","commonjs":"@gameworkers/react-simple-maps","commonjs2":"@gameworkers/react-simple-maps","amd":"@gameworkers/react-simple-maps"}
var react_simple_maps_ = __webpack_require__(2);

// EXTERNAL MODULE: external {"root":"NextReduxWrapper","commonjs":"next-redux-wrapper","commonjs2":"next-redux-wrapper","amd":"next-redux-wrapper"}
var external_root_NextReduxWrapper_commonjs_next_redux_wrapper_commonjs2_next_redux_wrapper_amd_next_redux_wrapper_ = __webpack_require__(5);
var external_root_NextReduxWrapper_commonjs_next_redux_wrapper_commonjs2_next_redux_wrapper_amd_next_redux_wrapper_default = /*#__PURE__*/__webpack_require__.n(external_root_NextReduxWrapper_commonjs_next_redux_wrapper_commonjs2_next_redux_wrapper_amd_next_redux_wrapper_);

// EXTERNAL MODULE: external {"root":"ReduxTooltip","commonjs":"redux-tooltip","commonjs2":"redux-tooltip","amd":"redux-tooltip"}
var external_root_ReduxTooltip_commonjs_redux_tooltip_commonjs2_redux_tooltip_amd_redux_tooltip_ = __webpack_require__(3);

// CONCATENATED MODULE: ./src/SvgContentElementWrapperWithDefs.js
function _objectWithoutPropertiesLoose(source, excluded) { if (source == null) return {}; var target = {}; var sourceKeys = Object.keys(source); var key, i; for (i = 0; i < sourceKeys.length; i++) { key = sourceKeys[i]; if (excluded.indexOf(key) >= 0) continue; target[key] = source[key]; } return target; }



var red = "#e5261a";

function MarkerPath(_ref) {
  var fill = _ref.fill,
      stroke = _ref.stroke,
      strokeWidth = _ref.strokeWidth;
  return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("path", {
    fill: fill,
    stroke: stroke,
    strokeWidth: strokeWidth,
    d: "\n        M-78.736,56.784c-89.734,0-162.725,73.003-162.725,162.732c0,36.571,25.572,94.673,76.01,172.686\n        c36.714,56.783,73.96,104.287,74.33,104.75c2.985,3.807,7.553,6.018,12.379,6.018c4.829,0,9.396-2.211,12.381-6.018\n        c0.37-0.463,37.622-47.971,74.33-104.75c50.438-78.013,76.012-136.121,76.012-172.686\n        C83.993,129.787,10.999,56.784-78.736,56.784z\n        M-140.173,241.81c-3.481,2.971-7.626,5.503-12.439,7.6\n        c-4.804,2.094-10.303,3.146-16.504,3.146c-5.293,0-10.114-0.871-14.47-2.616c-4.357-1.741-8.101-4.16-11.226-7.259\n        c-3.136-3.098-5.571-6.775-7.311-11.037c-1.742-4.259-2.616-8.903-2.616-13.936v-0.197c0-4.838,0.89-9.372,2.662-13.602\n        c1.772-4.222,4.248-7.937,7.405-11.131c3.16-3.193,6.909-5.712,11.234-7.55c4.317-1.839,9.061-2.758,14.228-2.758\n        c3.033,0,5.791,0.208,8.28,0.631c2.481,0.417,4.787,1.016,6.916,1.788c2.136,0.776,4.128,1.745,6.002,2.907\n        c1.876,1.158,3.682,2.481,5.421,3.966l-9.391,11.327c-1.291-1.098-2.583-2.067-3.868-2.906c-1.293-0.834-2.63-1.548-4.018-2.126\n        c-1.395-0.583-2.891-1.033-4.507-1.356c-1.614-0.322-3.387-0.484-5.316-0.484c-2.717,0-5.247,0.565-7.603,1.697\n        c-2.354,1.13-4.419,2.651-6.2,4.562c-1.771,1.906-3.158,4.139-4.158,6.692c-1.001,2.555-1.498,5.289-1.498,8.2v0.193\n        c0,3.104,0.497,5.968,1.498,8.588c1,2.617,2.418,4.898,4.263,6.838c1.835,1.941,4,3.44,6.481,4.513\n        c2.481,1.068,5.247,1.598,8.28,1.598c5.545,0,10.226-1.366,14.03-4.096v-9.749h-14.999v-12.872h29.423V241.81z\n        M-51.398,251.879\n        h-12.977l-15.29-44.337l-15.299,44.337h-12.974l-23.13-68.245h15.969l14.038,45.883l15.195-46.077h12.779l15.195,46.077\n        l14.039-45.883h15.58L-51.398,251.879z\n        M38.676,235.323c-1.395,3.81-3.381,6.988-5.956,9.536\n        c-2.585,2.553-5.712,4.455-9.39,5.711c-3.68,1.26-7.777,1.887-12.298,1.887c-9.099,0-16.262-2.517-21.49-7.551\n        c-5.224-5.033-7.839-12.55-7.839-22.554v-38.718h14.905v38.333c0,5.549,1.291,9.729,3.876,12.532\n        c2.575,2.809,6.161,4.217,10.744,4.217c4.579,0,8.161-1.359,10.746-4.07c2.575-2.71,3.868-6.775,3.868-12.194v-38.817h14.912\n        v38.234C40.755,227.032,40.061,231.518,38.676,235.323z\n      "
  });
}

function SvgContentElementWrapperWithDefs(_ref2) {
  var children = _ref2.children,
      forceGrayscale = _ref2.forceGrayscale,
      rest = _objectWithoutPropertiesLoose(_ref2, ["children", "forceGrayscale"]);

  return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.Fragment, null, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("defs", null, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("filter", {
    id: "grayscale"
  }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("feColorMatrix", {
    type: "saturate",
    values: "0"
  })), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("pattern", {
    id: "redpattern",
    patternUnits: "userSpaceOnUse",
    width: "8",
    height: "8"
  }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("rect", {
    style: {
      opacity: 0.7
    },
    width: "8",
    height: "8",
    fill: red
  })), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("pattern", {
    id: "lessredpattern",
    patternUnits: "userSpaceOnUse",
    width: "8",
    height: "8"
  }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("rect", {
    style: {
      opacity: 0.3
    },
    width: "8",
    height: "8",
    fill: red
  })), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("pattern", {
    id: "hardlyredpattern",
    patternUnits: "userSpaceOnUse",
    width: "8",
    height: "8"
  }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("rect", {
    style: {
      opacity: 0.04
    },
    width: "8",
    height: "8",
    fill: red
  })), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("g", {
    id: "mapmarker"
  }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("rect", {
    x: "-218.979",
    y: "163.179",
    fill: "#fff",
    width: "279.333",
    height: "106.667"
  }), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("g", null, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(MarkerPath, {
    stroke: "#fff",
    strokeWidth: 20
  }), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(MarkerPath, null)))), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("g", {
    filter: forceGrayscale ? "url(#grayscale)" : undefined
  }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.cloneElement(children, rest)));
}

SvgContentElementWrapperWithDefs.propTypes = {
  children: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.element.isRequired,
  forceGrayscale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.bool
};
/* harmony default export */ var src_SvgContentElementWrapperWithDefs = (SvgContentElementWrapperWithDefs);
// CONCATENATED MODULE: ./src/GWUMarker.js
function _extends() { _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; }; return _extends.apply(this, arguments); }

function GWUMarker_objectWithoutPropertiesLoose(source, excluded) { if (source == null) return {}; var target = {}; var sourceKeys = Object.keys(source); var key, i; for (i = 0; i < sourceKeys.length; i++) { key = sourceKeys[i]; if (excluded.indexOf(key) >= 0) continue; target[key] = source[key]; } return target; }

function _inheritsLoose(subClass, superClass) { subClass.prototype = Object.create(superClass.prototype); subClass.prototype.constructor = subClass; subClass.__proto__ = superClass; }





var GWUMarker_GWUMarker =
/*#__PURE__*/
function (_PureComponent) {
  _inheritsLoose(GWUMarker, _PureComponent);

  function GWUMarker() {
    return _PureComponent.apply(this, arguments) || this;
  }

  var _proto = GWUMarker.prototype;

  _proto.render = function render() {
    var _this$props = this.props,
        marker = _this$props.marker,
        scale = _this$props.scale,
        rest = GWUMarker_objectWithoutPropertiesLoose(_this$props, ["marker", "scale"]);

    return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(react_simple_maps_["Marker"], _extends({
      key: marker.name,
      marker: marker
    }, rest), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("g", {
      style: {
        cursor: "pointer"
      },
      transform: "scale(" + scale + "), translate(79, -528)"
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("use", {
      href: "#mapmarker"
    })));
  };

  return GWUMarker;
}(external_root_React_commonjs_react_commonjs2_react_amd_react_["PureComponent"]);

GWUMarker_GWUMarker.propTypes = {
  marker: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.shape({
    name: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.string.isRequired,
    // [lng, lat]
    coordinates: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.arrayOf(external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number).isRequired
  }).isRequired,
  scale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired
};
GWUMarker_GWUMarker.defaultProps = {
  scale: 0.09
};
/* harmony default export */ var src_GWUMarker = (GWUMarker_GWUMarker);
// EXTERNAL MODULE: external {"root":"Redux","commonjs":"redux","commonjs2":"redux","amd":"redux"}
var external_root_Redux_commonjs_redux_commonjs2_redux_amd_redux_ = __webpack_require__(4);

// CONCATENATED MODULE: ./src/tooltipStore.js


var initialState = {
  title: "With Redux Tooltip"
};

var appReducer = function appReducer(state, action) {
  if (state === void 0) {
    state = initialState;
  }

  switch (action.type) {
    default:
      return state;
  }
};

var tooltipStore_initStore = function initStore(initState) {
  if (initState === void 0) {
    initState = {
      appReducer: initialState
    };
  }

  return Object(external_root_Redux_commonjs_redux_commonjs2_redux_amd_redux_["createStore"])(Object(external_root_Redux_commonjs_redux_commonjs2_redux_amd_redux_["combineReducers"])({
    appReducer: appReducer,
    tooltip: external_root_ReduxTooltip_commonjs_redux_tooltip_commonjs2_redux_tooltip_amd_redux_tooltip_["reducer"]
  }), initState, window.__REDUX_DEVTOOLS_EXTENSION__ && window.__REDUX_DEVTOOLS_EXTENSION__());
};
// CONCATENATED MODULE: ./src/ChapterMapComponent.js
function _objectSpread(target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i] != null ? arguments[i] : {}; var ownKeys = Object.keys(source); if (typeof Object.getOwnPropertySymbols === 'function') { ownKeys = ownKeys.concat(Object.getOwnPropertySymbols(source).filter(function (sym) { return Object.getOwnPropertyDescriptor(source, sym).enumerable; })); } ownKeys.forEach(function (key) { _defineProperty(target, key, source[key]); }); } return target; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function ChapterMapComponent_inheritsLoose(subClass, superClass) { subClass.prototype = Object.create(superClass.prototype); subClass.prototype.constructor = subClass; subClass.__proto__ = superClass; }

function _assertThisInitialized(self) { if (self === void 0) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return self; }






external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.PropTypes = external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a; // for redux-tooltip compatibility

var wrapperStyles = {
  width: "100%",
  maxWidth: 980,
  margin: "0 auto",
  fontFamily: "Roboto, sans-serif"
};



var worldDataPathname = "https://gameworkers.github.io/data/third_party/world-50m.json";
var smallWorldDataPathname = "https://gameworkers.github.io/data/third_party/world-110m.json";
var membersPathname = "https://gameworkers.github.io/data/members.json";
var worldDataPromise;
var smallWorldDataPromise;
var membersPromise;

function getWorldData() {
  return worldDataPromise = worldDataPromise || fetch(worldDataPathname).then(function (res) {
    return res.json();
  });
}

function getSmallWorldData() {
  return smallWorldDataPromise = smallWorldDataPromise || fetch(smallWorldDataPathname).then(function (res) {
    return res.json();
  });
}

function getMembers() {
  return membersPromise = membersPromise || fetch(membersPathname).then(function (res) {
    return res.json();
  });
}

var countryNameKeys = ["ABBREV", "CONTINENT", "FORMAL_EN", "ISO_A2", "ISO_A3", "NAME", "NAME_LONG"];

function geographyMatchesCountryString(geography, countryString) {
  return countryNameKeys.some(function (key) {
    return geography.properties[key] === countryString;
  });
}

var ChapterMapComponent_ChapterMapComponent =
/*#__PURE__*/
function (_PureComponent) {
  ChapterMapComponent_inheritsLoose(ChapterMapComponent, _PureComponent);

  function ChapterMapComponent(props) {
    var _this;

    _this = _PureComponent.call(this, props) || this;
    _this.state = {
      worldData: null,
      members: null,
      markers: null
    };
    _this.focusedMarker = null;
    _this.handleCountryMove = _this.handleCountryMove.bind(_assertThisInitialized(_assertThisInitialized(_this)));
    _this.handleCountryLeave = _this.handleCountryLeave.bind(_assertThisInitialized(_assertThisInitialized(_this)));
    _this.handleMarkerMove = _this.handleMarkerMove.bind(_assertThisInitialized(_assertThisInitialized(_this)));
    _this.handleMarkerLeave = _this.handleMarkerLeave.bind(_assertThisInitialized(_assertThisInitialized(_this)));
    _this.handleMarkerClick = _this.handleMarkerClick.bind(_assertThisInitialized(_assertThisInitialized(_this)));
    _this.handleOutsideMarkerClick = _this.handleOutsideMarkerClick.bind(_assertThisInitialized(_assertThisInitialized(_this)));
    return _this;
  }

  var _proto = ChapterMapComponent.prototype;

  _proto.componentDidMount = function componentDidMount() {
    var _this2 = this;

    getWorldData().then(function (worldData) {
      _this2.setState({
        worldData: worldData
      });
    });
    getSmallWorldData().then(function (worldData) {
      _this2.setState(function (state) {
        return state.worldData ? null : {
          worldData: worldData
        };
      });
    });
    getMembers().then(function (members) {
      _this2.setState({
        members: members,
        markers: members.filter(function (member) {
          return member.isChapter;
        }).sort(function (a, b) {
          // we want to make sure markers lower on the map are painted in front
          if (a.lat < b.lat) {
            return 1;
          }

          if (a.lat > b.lat) {
            return -1;
          }

          if (a.lng < b.lng) {
            return 1;
          }

          if (a.lng > b.lng) {
            return -1;
          }

          return 0;
        }).map(function (member) {
          return {
            name: member.location,
            coordinates: [member.lng, member.lat],
            data: member
          };
        })
      });
    });
    window.addEventListener("click", this.handleOutsideMarkerClick);
  };

  _proto.componentWillUnmount = function componentWillUnmount() {
    window.removeEventListener("click", this.handleOutsideMarkerClick);
  };

  _proto.handleCountryMove = function handleCountryMove(geography, evt) {
    if (this.focusedMarker) {
      return;
    }

    this.dispatchTooltip(evt, geography.properties.NAME_LONG);
  };

  _proto.handleCountryLeave = function handleCountryLeave() {
    if (this.focusedMarker) {
      return;
    }

    this.hideTooltip();
  };

  _proto.handleMarkerMove = function handleMarkerMove(marker, evt) {
    if (this.focusedMarker) {
      return;
    }

    this.dispatchTooltip(evt, marker.name);
  };

  _proto.handleMarkerLeave = function handleMarkerLeave() {
    if (this.focusedMarker) {
      return;
    }

    this.hideTooltip();
  };

  _proto.handleMarkerClick = function handleMarkerClick(marker, coords, evt) {
    evt.stopPropagation();
    evt.nativeEvent.stopImmediatePropagation();
    this.focusedMarker = this.focusedMarker === marker ? null : marker;

    if (this.focusedMarker) {
      var _marker$data = marker.data,
          location = _marker$data.location,
          chapterInfo = _marker$data.chapterInfo;
      var content = location;

      if (chapterInfo) {
        var description = chapterInfo.description,
            applicationLink = chapterInfo.applicationLink,
            twitter = chapterInfo.twitter,
            email = chapterInfo.email,
            website = chapterInfo.website;
        content = "<h3>" + location + "</h3>";

        if (description) {
          content += "<p>" + description + "</p>";
        }

        if (twitter || email || website) {
          content += "<p>";
        }

        if (twitter) {
          content += "\n            <strong>Twitter:</strong>\n            <a href=\"https://twitter.com/" + twitter + "\">\n              @" + twitter + "\n            </a>\n            <br />\n          ";
        }

        if (email) {
          content += "\n            <strong>Email:</strong>\n            <a href=\"mailto:" + email + "\">\n              " + email + "\n            </a>\n            <br />\n          ";
        }

        if (website) {
          content += "\n            <strong>Website:</strong>\n            <a href=\"" + (website.indexOf("http") === 0 ? "" : "http://") + website + "\">\n              " + website + "\n            </a>\n            <br />\n          ";
        }

        if (twitter || email || website) {
          content += "</p>";
        }

        if (applicationLink) {
          content += "\n            <p>\n              <a href=\"" + applicationLink + "\">\n                Apply here!\n              </a>\n            </p>\n          ";
        }

        content = "<div class=\"" + this.props.tooltipClassName + "\">" + content + "</div>";
      }

      this.dispatchTooltip(evt, content);
    } else {
      this.hideTooltip();
    }
  };

  _proto.handleOutsideMarkerClick = function handleOutsideMarkerClick() {
    this.hideTooltip();
  };

  _proto.dispatchTooltip = function dispatchTooltip(evt, content) {
    var x = evt.clientX;
    var y = evt.clientY + window.pageYOffset;
    this.props.dispatch(external_root_ReduxTooltip_commonjs_redux_tooltip_commonjs2_redux_tooltip_amd_redux_tooltip_["actions"].show({
      origin: {
        x: x,
        y: y
      },
      content: content
    }));
  };

  _proto.hideTooltip = function hideTooltip() {
    var _this3 = this;

    this.props.dispatch(external_root_ReduxTooltip_commonjs_redux_tooltip_commonjs2_redux_tooltip_amd_redux_tooltip_["actions"].hide());
    Promise.resolve().then(function () {
      _this3.focusedMarker = null;
    });
  };

  _proto.render = function render() {
    var _this4 = this;

    var _this$props = this.props,
        centerLat = _this$props.centerLat,
        centerLng = _this$props.centerLng,
        width = _this$props.width,
        height = _this$props.height,
        scale = _this$props.scale,
        isGeographyIncluded = _this$props.isGeographyIncluded,
        markerScale = _this$props.markerScale,
        style = _this$props.style,
        className = _this$props.className,
        forceGrayscale = _this$props.forceGrayscale,
        zoom = _this$props.zoom,
        enablePanning = _this$props.enablePanning,
        projection = _this$props.projection;
    var _this$state = this.state,
        worldData = _this$state.worldData,
        members = _this$state.members,
        markers = _this$state.markers;
    var loading = !(worldData && members && markers);
    return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("div", {
      className: className,
      style: _objectSpread({}, style || {}, {
        width: width,
        height: height
      })
    }, loading && external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("div", null, "Loading..."), !loading && external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.Fragment, null, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(react_simple_maps_["ComposableMap"], {
      projection: projection,
      projectionConfig: {
        scale: scale
      },
      width: width,
      height: height,
      style: {
        width: "100%",
        height: "auto",
        backgroundColor: "#fff"
      }
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(src_SvgContentElementWrapperWithDefs, {
      forceGrayscale: forceGrayscale
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(react_simple_maps_["ZoomableGroup"], {
      center: [centerLng, centerLat],
      disablePanning: !enablePanning,
      zoom: zoom
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(react_simple_maps_["Geographies"], {
      geography: this.state.worldData
    }, function (geographies, projection) {
      return geographies.filter(isGeographyIncluded).map(function (geography, i) {
        var hasMatchingPoint = members.some(function (member) {
          return geographyMatchesCountryString(geography, member.country);
        });
        var style = {
          fill: hasMatchingPoint ? "url(#redpattern)" : "url(#hardlyredpattern)",
          stroke: "#222",
          strokeWidth: 0.5 / zoom,
          outline: "none"
        };
        return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(react_simple_maps_["Geography"], {
          key: i,
          geography: geography,
          projection: projection,
          style: {
            default: style,
            hover: style,
            pressed: style
          },
          onMouseMove: hasMatchingPoint && _this4.handleCountryMove,
          onMouseLeave: hasMatchingPoint && _this4.handleCountryLeave
        });
      });
    }), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(react_simple_maps_["Markers"], null, markers.map(function (marker) {
      return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(src_GWUMarker, {
        key: marker.name,
        marker: marker,
        scale: markerScale,
        onClick: _this4.handleMarkerClick,
        onMouseMove: _this4.handleMarkerMove,
        onMouseLeave: _this4.handleMarkerLeave
      });
    }))))), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(external_root_ReduxTooltip_commonjs_redux_tooltip_commonjs2_redux_tooltip_amd_redux_tooltip_["Tooltip"], null)));
  };

  return ChapterMapComponent;
}(external_root_React_commonjs_react_commonjs2_react_amd_react_["PureComponent"]);

ChapterMapComponent_ChapterMapComponent = external_root_NextReduxWrapper_commonjs_next_redux_wrapper_commonjs2_next_redux_wrapper_amd_next_redux_wrapper_default()(tooltipStore_initStore)(ChapterMapComponent_ChapterMapComponent);
ChapterMapComponent_ChapterMapComponent.displayName = "ChapterMapComponent";
ChapterMapComponent_ChapterMapComponent.propTypes = {
  centerLat: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  centerLng: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  width: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  height: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  scale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  isGeographyIncluded: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.func.isRequired,
  markerScale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  forceGrayscale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.bool.isRequired,
  tooltipClassName: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.string.isRequired,
  zoom: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  enablePanning: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.bool.isRequired,
  projection: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.string
};
ChapterMapComponent_ChapterMapComponent.defaultProps = {
  centerLat: 0,
  centerLng: 0,
  width: 980,
  height: 551,
  scale: 205,
  isGeographyIncluded: function isGeographyIncluded() {
    return true;
  },
  markerScale: 0.09,
  forceGrayscale: false,
  tooltipClassName: "gwu_chapter_tooltip",
  zoom: 1,
  enablePanning: false,
  projection: "times"
};
/* harmony default export */ var src_ChapterMapComponent = (ChapterMapComponent_ChapterMapComponent);
// CONCATENATED MODULE: ./src/index.js

/* harmony default export */ var src = __webpack_exports__["default"] = (src_ChapterMapComponent);

/***/ })
/******/ ])["default"];
});
//# sourceMappingURL=chapter-map-component.js.map