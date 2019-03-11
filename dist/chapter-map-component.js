(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory(require("react"), require("prop-types"));
	else if(typeof define === 'function' && define.amd)
		define(["react", "prop-types"], factory);
	else if(typeof exports === 'object')
		exports["ChapterMapComponent"] = factory(require("react"), require("prop-types"));
	else
		root["ChapterMapComponent"] = factory(root["React"], root["PropTypes"]);
})((typeof self !== "undefined" ? self : this), function(__WEBPACK_EXTERNAL_MODULE__1__, __WEBPACK_EXTERNAL_MODULE__3__) {
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
/******/ 	return __webpack_require__(__webpack_require__.s = 23);
/******/ })
/************************************************************************/
/******/ ([
/* 0 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
__webpack_require__.r(__webpack_exports__);

// CONCATENATED MODULE: ./node_modules/d3-geo/src/adder.js
// Adds floating point numbers with twice the normal precision.
// Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and
// Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)
// 305–363 (1997).
// Code adapted from GeographicLib by Charles F. F. Karney,
// http://geographiclib.sourceforge.net/

/* harmony default export */ var adder = (function() {
  return new Adder;
});

function Adder() {
  this.reset();
}

Adder.prototype = {
  constructor: Adder,
  reset: function() {
    this.s = // rounded value
    this.t = 0; // exact error
  },
  add: function(y) {
    add(temp, y, this.t);
    add(this, temp.s, this.s);
    if (this.s) this.t += temp.t;
    else this.s = temp.t;
  },
  valueOf: function() {
    return this.s;
  }
};

var temp = new Adder;

function add(adder, a, b) {
  var x = adder.s = a + b,
      bv = x - a,
      av = x - bv;
  adder.t = (a - av) + (b - bv);
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/math.js
var epsilon = 1e-6;
var epsilon2 = 1e-12;
var pi = Math.PI;
var halfPi = pi / 2;
var quarterPi = pi / 4;
var tau = pi * 2;

var degrees = 180 / pi;
var radians = pi / 180;

var abs = Math.abs;
var atan = Math.atan;
var atan2 = Math.atan2;
var cos = Math.cos;
var ceil = Math.ceil;
var exp = Math.exp;
var floor = Math.floor;
var log = Math.log;
var pow = Math.pow;
var sin = Math.sin;
var math_sign = Math.sign || function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; };
var sqrt = Math.sqrt;
var tan = Math.tan;

function acos(x) {
  return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
}

function asin(x) {
  return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
}

function haversin(x) {
  return (x = sin(x / 2)) * x;
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/noop.js
function noop() {}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/stream.js
function streamGeometry(geometry, stream) {
  if (geometry && streamGeometryType.hasOwnProperty(geometry.type)) {
    streamGeometryType[geometry.type](geometry, stream);
  }
}

var streamObjectType = {
  Feature: function(object, stream) {
    streamGeometry(object.geometry, stream);
  },
  FeatureCollection: function(object, stream) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n) streamGeometry(features[i].geometry, stream);
  }
};

var streamGeometryType = {
  Sphere: function(object, stream) {
    stream.sphere();
  },
  Point: function(object, stream) {
    object = object.coordinates;
    stream.point(object[0], object[1], object[2]);
  },
  MultiPoint: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) object = coordinates[i], stream.point(object[0], object[1], object[2]);
  },
  LineString: function(object, stream) {
    streamLine(object.coordinates, stream, 0);
  },
  MultiLineString: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) streamLine(coordinates[i], stream, 0);
  },
  Polygon: function(object, stream) {
    streamPolygon(object.coordinates, stream);
  },
  MultiPolygon: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) streamPolygon(coordinates[i], stream);
  },
  GeometryCollection: function(object, stream) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n) streamGeometry(geometries[i], stream);
  }
};

function streamLine(coordinates, stream, closed) {
  var i = -1, n = coordinates.length - closed, coordinate;
  stream.lineStart();
  while (++i < n) coordinate = coordinates[i], stream.point(coordinate[0], coordinate[1], coordinate[2]);
  stream.lineEnd();
}

function streamPolygon(coordinates, stream) {
  var i = -1, n = coordinates.length;
  stream.polygonStart();
  while (++i < n) streamLine(coordinates[i], stream, 1);
  stream.polygonEnd();
}

/* harmony default export */ var src_stream = (function(object, stream) {
  if (object && streamObjectType.hasOwnProperty(object.type)) {
    streamObjectType[object.type](object, stream);
  } else {
    streamGeometry(object, stream);
  }
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/area.js





var areaRingSum = adder();

var areaSum = adder(),
    area_lambda00,
    phi00,
    area_lambda0,
    area_cosPhi0,
    area_sinPhi0;

var areaStream = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    areaRingSum.reset();
    areaStream.lineStart = areaRingStart;
    areaStream.lineEnd = areaRingEnd;
  },
  polygonEnd: function() {
    var areaRing = +areaRingSum;
    areaSum.add(areaRing < 0 ? tau + areaRing : areaRing);
    this.lineStart = this.lineEnd = this.point = noop;
  },
  sphere: function() {
    areaSum.add(tau);
  }
};

function areaRingStart() {
  areaStream.point = areaPointFirst;
}

function areaRingEnd() {
  areaPoint(area_lambda00, phi00);
}

function areaPointFirst(lambda, phi) {
  areaStream.point = areaPoint;
  area_lambda00 = lambda, phi00 = phi;
  lambda *= radians, phi *= radians;
  area_lambda0 = lambda, area_cosPhi0 = cos(phi = phi / 2 + quarterPi), area_sinPhi0 = sin(phi);
}

function areaPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  phi = phi / 2 + quarterPi; // half the angular distance from south pole

  // Spherical excess E for a spherical triangle with vertices: south pole,
  // previous point, current point.  Uses a formula derived from Cagnoli’s
  // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).
  var dLambda = lambda - area_lambda0,
      sdLambda = dLambda >= 0 ? 1 : -1,
      adLambda = sdLambda * dLambda,
      cosPhi = cos(phi),
      sinPhi = sin(phi),
      k = area_sinPhi0 * sinPhi,
      u = area_cosPhi0 * cosPhi + k * cos(adLambda),
      v = k * sdLambda * sin(adLambda);
  areaRingSum.add(atan2(v, u));

  // Advance the previous points.
  area_lambda0 = lambda, area_cosPhi0 = cosPhi, area_sinPhi0 = sinPhi;
}

/* harmony default export */ var src_area = (function(object) {
  areaSum.reset();
  src_stream(object, areaStream);
  return areaSum * 2;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/cartesian.js


function cartesian_spherical(cartesian) {
  return [atan2(cartesian[1], cartesian[0]), asin(cartesian[2])];
}

function cartesian_cartesian(spherical) {
  var lambda = spherical[0], phi = spherical[1], cosPhi = cos(phi);
  return [cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi)];
}

function cartesianDot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cartesianCross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

// TODO return a
function cartesianAddInPlace(a, b) {
  a[0] += b[0], a[1] += b[1], a[2] += b[2];
}

function cartesianScale(vector, k) {
  return [vector[0] * k, vector[1] * k, vector[2] * k];
}

// TODO return d
function cartesianNormalizeInPlace(d) {
  var l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  d[0] /= l, d[1] /= l, d[2] /= l;
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/bounds.js






var bounds_lambda0, bounds_phi0, bounds_lambda1, bounds_phi1, // bounds
    bounds_lambda2, // previous lambda-coordinate
    bounds_lambda00, bounds_phi00, // first point
    bounds_p0, // previous 3D point
    deltaSum = adder(),
    ranges,
    range;

var boundsStream = {
  point: boundsPoint,
  lineStart: boundsLineStart,
  lineEnd: boundsLineEnd,
  polygonStart: function() {
    boundsStream.point = boundsRingPoint;
    boundsStream.lineStart = boundsRingStart;
    boundsStream.lineEnd = boundsRingEnd;
    deltaSum.reset();
    areaStream.polygonStart();
  },
  polygonEnd: function() {
    areaStream.polygonEnd();
    boundsStream.point = boundsPoint;
    boundsStream.lineStart = boundsLineStart;
    boundsStream.lineEnd = boundsLineEnd;
    if (areaRingSum < 0) bounds_lambda0 = -(bounds_lambda1 = 180), bounds_phi0 = -(bounds_phi1 = 90);
    else if (deltaSum > epsilon) bounds_phi1 = 90;
    else if (deltaSum < -epsilon) bounds_phi0 = -90;
    range[0] = bounds_lambda0, range[1] = bounds_lambda1;
  }
};

function boundsPoint(lambda, phi) {
  ranges.push(range = [bounds_lambda0 = lambda, bounds_lambda1 = lambda]);
  if (phi < bounds_phi0) bounds_phi0 = phi;
  if (phi > bounds_phi1) bounds_phi1 = phi;
}

function bounds_linePoint(lambda, phi) {
  var p = cartesian_cartesian([lambda * radians, phi * radians]);
  if (bounds_p0) {
    var normal = cartesianCross(bounds_p0, p),
        equatorial = [normal[1], -normal[0], 0],
        inflection = cartesianCross(equatorial, normal);
    cartesianNormalizeInPlace(inflection);
    inflection = cartesian_spherical(inflection);
    var delta = lambda - bounds_lambda2,
        sign = delta > 0 ? 1 : -1,
        lambdai = inflection[0] * degrees * sign,
        phii,
        antimeridian = abs(delta) > 180;
    if (antimeridian ^ (sign * bounds_lambda2 < lambdai && lambdai < sign * lambda)) {
      phii = inflection[1] * degrees;
      if (phii > bounds_phi1) bounds_phi1 = phii;
    } else if (lambdai = (lambdai + 360) % 360 - 180, antimeridian ^ (sign * bounds_lambda2 < lambdai && lambdai < sign * lambda)) {
      phii = -inflection[1] * degrees;
      if (phii < bounds_phi0) bounds_phi0 = phii;
    } else {
      if (phi < bounds_phi0) bounds_phi0 = phi;
      if (phi > bounds_phi1) bounds_phi1 = phi;
    }
    if (antimeridian) {
      if (lambda < bounds_lambda2) {
        if (bounds_angle(bounds_lambda0, lambda) > bounds_angle(bounds_lambda0, bounds_lambda1)) bounds_lambda1 = lambda;
      } else {
        if (bounds_angle(lambda, bounds_lambda1) > bounds_angle(bounds_lambda0, bounds_lambda1)) bounds_lambda0 = lambda;
      }
    } else {
      if (bounds_lambda1 >= bounds_lambda0) {
        if (lambda < bounds_lambda0) bounds_lambda0 = lambda;
        if (lambda > bounds_lambda1) bounds_lambda1 = lambda;
      } else {
        if (lambda > bounds_lambda2) {
          if (bounds_angle(bounds_lambda0, lambda) > bounds_angle(bounds_lambda0, bounds_lambda1)) bounds_lambda1 = lambda;
        } else {
          if (bounds_angle(lambda, bounds_lambda1) > bounds_angle(bounds_lambda0, bounds_lambda1)) bounds_lambda0 = lambda;
        }
      }
    }
  } else {
    ranges.push(range = [bounds_lambda0 = lambda, bounds_lambda1 = lambda]);
  }
  if (phi < bounds_phi0) bounds_phi0 = phi;
  if (phi > bounds_phi1) bounds_phi1 = phi;
  bounds_p0 = p, bounds_lambda2 = lambda;
}

function boundsLineStart() {
  boundsStream.point = bounds_linePoint;
}

function boundsLineEnd() {
  range[0] = bounds_lambda0, range[1] = bounds_lambda1;
  boundsStream.point = boundsPoint;
  bounds_p0 = null;
}

function boundsRingPoint(lambda, phi) {
  if (bounds_p0) {
    var delta = lambda - bounds_lambda2;
    deltaSum.add(abs(delta) > 180 ? delta + (delta > 0 ? 360 : -360) : delta);
  } else {
    bounds_lambda00 = lambda, bounds_phi00 = phi;
  }
  areaStream.point(lambda, phi);
  bounds_linePoint(lambda, phi);
}

function boundsRingStart() {
  areaStream.lineStart();
}

function boundsRingEnd() {
  boundsRingPoint(bounds_lambda00, bounds_phi00);
  areaStream.lineEnd();
  if (abs(deltaSum) > epsilon) bounds_lambda0 = -(bounds_lambda1 = 180);
  range[0] = bounds_lambda0, range[1] = bounds_lambda1;
  bounds_p0 = null;
}

// Finds the left-right distance between two longitudes.
// This is almost the same as (lambda1 - lambda0 + 360°) % 360°, except that we want
// the distance between ±180° to be 360°.
function bounds_angle(lambda0, lambda1) {
  return (lambda1 -= lambda0) < 0 ? lambda1 + 360 : lambda1;
}

function rangeCompare(a, b) {
  return a[0] - b[0];
}

function rangeContains(range, x) {
  return range[0] <= range[1] ? range[0] <= x && x <= range[1] : x < range[0] || range[1] < x;
}

/* harmony default export */ var bounds = (function(feature) {
  var i, n, a, b, merged, deltaMax, delta;

  bounds_phi1 = bounds_lambda1 = -(bounds_lambda0 = bounds_phi0 = Infinity);
  ranges = [];
  src_stream(feature, boundsStream);

  // First, sort ranges by their minimum longitudes.
  if (n = ranges.length) {
    ranges.sort(rangeCompare);

    // Then, merge any ranges that overlap.
    for (i = 1, a = ranges[0], merged = [a]; i < n; ++i) {
      b = ranges[i];
      if (rangeContains(a, b[0]) || rangeContains(a, b[1])) {
        if (bounds_angle(a[0], b[1]) > bounds_angle(a[0], a[1])) a[1] = b[1];
        if (bounds_angle(b[0], a[1]) > bounds_angle(a[0], a[1])) a[0] = b[0];
      } else {
        merged.push(a = b);
      }
    }

    // Finally, find the largest gap between the merged ranges.
    // The final bounding box will be the inverse of this gap.
    for (deltaMax = -Infinity, n = merged.length - 1, i = 0, a = merged[n]; i <= n; a = b, ++i) {
      b = merged[i];
      if ((delta = bounds_angle(a[1], b[0])) > deltaMax) deltaMax = delta, bounds_lambda0 = b[0], bounds_lambda1 = a[1];
    }
  }

  ranges = range = null;

  return bounds_lambda0 === Infinity || bounds_phi0 === Infinity
      ? [[NaN, NaN], [NaN, NaN]]
      : [[bounds_lambda0, bounds_phi0], [bounds_lambda1, bounds_phi1]];
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/centroid.js




var W0, W1,
    centroid_X0, centroid_Y0, Z0,
    centroid_X1, centroid_Y1, Z1,
    X2, Y2, Z2,
    centroid_lambda00, centroid_phi00, // first point
    centroid_x0, centroid_y0, z0; // previous point

var centroidStream = {
  sphere: noop,
  point: centroidPoint,
  lineStart: centroidLineStart,
  lineEnd: centroidLineEnd,
  polygonStart: function() {
    centroidStream.lineStart = centroidRingStart;
    centroidStream.lineEnd = centroidRingEnd;
  },
  polygonEnd: function() {
    centroidStream.lineStart = centroidLineStart;
    centroidStream.lineEnd = centroidLineEnd;
  }
};

// Arithmetic mean of Cartesian vectors.
function centroidPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  centroidPointCartesian(cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi));
}

function centroidPointCartesian(x, y, z) {
  ++W0;
  centroid_X0 += (x - centroid_X0) / W0;
  centroid_Y0 += (y - centroid_Y0) / W0;
  Z0 += (z - Z0) / W0;
}

function centroidLineStart() {
  centroidStream.point = centroidLinePointFirst;
}

function centroidLinePointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  centroid_x0 = cosPhi * cos(lambda);
  centroid_y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidStream.point = centroidLinePoint;
  centroidPointCartesian(centroid_x0, centroid_y0, z0);
}

function centroidLinePoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi),
      x = cosPhi * cos(lambda),
      y = cosPhi * sin(lambda),
      z = sin(phi),
      w = atan2(sqrt((w = centroid_y0 * z - z0 * y) * w + (w = z0 * x - centroid_x0 * z) * w + (w = centroid_x0 * y - centroid_y0 * x) * w), centroid_x0 * x + centroid_y0 * y + z0 * z);
  W1 += w;
  centroid_X1 += w * (centroid_x0 + (centroid_x0 = x));
  centroid_Y1 += w * (centroid_y0 + (centroid_y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(centroid_x0, centroid_y0, z0);
}

function centroidLineEnd() {
  centroidStream.point = centroidPoint;
}

// See J. E. Brock, The Inertia Tensor for a Spherical Triangle,
// J. Applied Mechanics 42, 239 (1975).
function centroidRingStart() {
  centroidStream.point = centroidRingPointFirst;
}

function centroidRingEnd() {
  centroidRingPoint(centroid_lambda00, centroid_phi00);
  centroidStream.point = centroidPoint;
}

function centroidRingPointFirst(lambda, phi) {
  centroid_lambda00 = lambda, centroid_phi00 = phi;
  lambda *= radians, phi *= radians;
  centroidStream.point = centroidRingPoint;
  var cosPhi = cos(phi);
  centroid_x0 = cosPhi * cos(lambda);
  centroid_y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidPointCartesian(centroid_x0, centroid_y0, z0);
}

function centroidRingPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi),
      x = cosPhi * cos(lambda),
      y = cosPhi * sin(lambda),
      z = sin(phi),
      cx = centroid_y0 * z - z0 * y,
      cy = z0 * x - centroid_x0 * z,
      cz = centroid_x0 * y - centroid_y0 * x,
      m = sqrt(cx * cx + cy * cy + cz * cz),
      w = asin(m), // line weight = angle
      v = m && -w / m; // area weight multiplier
  X2 += v * cx;
  Y2 += v * cy;
  Z2 += v * cz;
  W1 += w;
  centroid_X1 += w * (centroid_x0 + (centroid_x0 = x));
  centroid_Y1 += w * (centroid_y0 + (centroid_y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(centroid_x0, centroid_y0, z0);
}

/* harmony default export */ var centroid = (function(object) {
  W0 = W1 =
  centroid_X0 = centroid_Y0 = Z0 =
  centroid_X1 = centroid_Y1 = Z1 =
  X2 = Y2 = Z2 = 0;
  src_stream(object, centroidStream);

  var x = X2,
      y = Y2,
      z = Z2,
      m = x * x + y * y + z * z;

  // If the area-weighted ccentroid is undefined, fall back to length-weighted ccentroid.
  if (m < epsilon2) {
    x = centroid_X1, y = centroid_Y1, z = Z1;
    // If the feature has zero length, fall back to arithmetic mean of point vectors.
    if (W1 < epsilon) x = centroid_X0, y = centroid_Y0, z = Z0;
    m = x * x + y * y + z * z;
    // If the feature still has an undefined ccentroid, then return.
    if (m < epsilon2) return [NaN, NaN];
  }

  return [atan2(y, x) * degrees, asin(z / sqrt(m)) * degrees];
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/constant.js
/* harmony default export */ var constant = (function(x) {
  return function() {
    return x;
  };
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/compose.js
/* harmony default export */ var compose = (function(a, b) {

  function compose(x, y) {
    return x = a(x, y), b(x[0], x[1]);
  }

  if (a.invert && b.invert) compose.invert = function(x, y) {
    return x = b.invert(x, y), x && a.invert(x[0], x[1]);
  };

  return compose;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/rotation.js



function rotationIdentity(lambda, phi) {
  return [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
}

rotationIdentity.invert = rotationIdentity;

function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
  return (deltaLambda %= tau) ? (deltaPhi || deltaGamma ? compose(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma))
    : rotationLambda(deltaLambda))
    : (deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma)
    : rotationIdentity);
}

function forwardRotationLambda(deltaLambda) {
  return function(lambda, phi) {
    return lambda += deltaLambda, [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
  };
}

function rotationLambda(deltaLambda) {
  var rotation = forwardRotationLambda(deltaLambda);
  rotation.invert = forwardRotationLambda(-deltaLambda);
  return rotation;
}

function rotationPhiGamma(deltaPhi, deltaGamma) {
  var cosDeltaPhi = cos(deltaPhi),
      sinDeltaPhi = sin(deltaPhi),
      cosDeltaGamma = cos(deltaGamma),
      sinDeltaGamma = sin(deltaGamma);

  function rotation(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi),
        k = z * cosDeltaPhi + x * sinDeltaPhi;
    return [
      atan2(y * cosDeltaGamma - k * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi),
      asin(k * cosDeltaGamma + y * sinDeltaGamma)
    ];
  }

  rotation.invert = function(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi),
        k = z * cosDeltaGamma - y * sinDeltaGamma;
    return [
      atan2(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k * sinDeltaPhi),
      asin(k * cosDeltaPhi - x * sinDeltaPhi)
    ];
  };

  return rotation;
}

/* harmony default export */ var src_rotation = (function(rotate) {
  rotate = rotateRadians(rotate[0] * radians, rotate[1] * radians, rotate.length > 2 ? rotate[2] * radians : 0);

  function forward(coordinates) {
    coordinates = rotate(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  }

  forward.invert = function(coordinates) {
    coordinates = rotate.invert(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  };

  return forward;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/circle.js





// Generates a circle centered at [0°, 0°], with a given radius and precision.
function circleStream(stream, radius, delta, direction, t0, t1) {
  if (!delta) return;
  var cosRadius = cos(radius),
      sinRadius = sin(radius),
      step = direction * delta;
  if (t0 == null) {
    t0 = radius + direction * tau;
    t1 = radius - step / 2;
  } else {
    t0 = circleRadius(cosRadius, t0);
    t1 = circleRadius(cosRadius, t1);
    if (direction > 0 ? t0 < t1 : t0 > t1) t0 += direction * tau;
  }
  for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
    point = cartesian_spherical([cosRadius, -sinRadius * cos(t), -sinRadius * sin(t)]);
    stream.point(point[0], point[1]);
  }
}

// Returns the signed angle of a cartesian point relative to [cosRadius, 0, 0].
function circleRadius(cosRadius, point) {
  point = cartesian_cartesian(point), point[0] -= cosRadius;
  cartesianNormalizeInPlace(point);
  var radius = acos(-point[1]);
  return ((-point[2] < 0 ? -radius : radius) + tau - epsilon) % tau;
}

/* harmony default export */ var src_circle = (function() {
  var center = constant([0, 0]),
      radius = constant(90),
      precision = constant(6),
      ring,
      rotate,
      stream = {point: point};

  function point(x, y) {
    ring.push(x = rotate(x, y));
    x[0] *= degrees, x[1] *= degrees;
  }

  function circle() {
    var c = center.apply(this, arguments),
        r = radius.apply(this, arguments) * radians,
        p = precision.apply(this, arguments) * radians;
    ring = [];
    rotate = rotateRadians(-c[0] * radians, -c[1] * radians, 0).invert;
    circleStream(stream, r, p, 1);
    c = {type: "Polygon", coordinates: [ring]};
    ring = rotate = null;
    return c;
  }

  circle.center = function(_) {
    return arguments.length ? (center = typeof _ === "function" ? _ : constant([+_[0], +_[1]]), circle) : center;
  };

  circle.radius = function(_) {
    return arguments.length ? (radius = typeof _ === "function" ? _ : constant(+_), circle) : radius;
  };

  circle.precision = function(_) {
    return arguments.length ? (precision = typeof _ === "function" ? _ : constant(+_), circle) : precision;
  };

  return circle;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/buffer.js


/* harmony default export */ var buffer = (function() {
  var lines = [],
      line;
  return {
    point: function(x, y) {
      line.push([x, y]);
    },
    lineStart: function() {
      lines.push(line = []);
    },
    lineEnd: noop,
    rejoin: function() {
      if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));
    },
    result: function() {
      var result = lines;
      lines = [];
      line = null;
      return result;
    }
  };
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/line.js
/* harmony default export */ var clip_line = (function(a, b, x0, y0, x1, y1) {
  var ax = a[0],
      ay = a[1],
      bx = b[0],
      by = b[1],
      t0 = 0,
      t1 = 1,
      dx = bx - ax,
      dy = by - ay,
      r;

  r = x0 - ax;
  if (!dx && r > 0) return;
  r /= dx;
  if (dx < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dx > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = x1 - ax;
  if (!dx && r < 0) return;
  r /= dx;
  if (dx < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dx > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  r = y0 - ay;
  if (!dy && r > 0) return;
  r /= dy;
  if (dy < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dy > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = y1 - ay;
  if (!dy && r < 0) return;
  r /= dy;
  if (dy < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dy > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  if (t0 > 0) a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
  if (t1 < 1) b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
  return true;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/pointEqual.js


/* harmony default export */ var pointEqual = (function(a, b) {
  return abs(a[0] - b[0]) < epsilon && abs(a[1] - b[1]) < epsilon;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/polygon.js


function Intersection(point, points, other, entry) {
  this.x = point;
  this.z = points;
  this.o = other; // another intersection
  this.e = entry; // is an entry?
  this.v = false; // visited
  this.n = this.p = null; // next & previous
}

// A generalized polygon clipping algorithm: given a polygon that has been cut
// into its visible line segments, and rejoins the segments by interpolating
// along the clip edge.
/* harmony default export */ var clip_polygon = (function(segments, compareIntersection, startInside, interpolate, stream) {
  var subject = [],
      clip = [],
      i,
      n;

  segments.forEach(function(segment) {
    if ((n = segment.length - 1) <= 0) return;
    var n, p0 = segment[0], p1 = segment[n], x;

    // If the first and last points of a segment are coincident, then treat as a
    // closed ring. TODO if all rings are closed, then the winding order of the
    // exterior ring should be checked.
    if (pointEqual(p0, p1)) {
      stream.lineStart();
      for (i = 0; i < n; ++i) stream.point((p0 = segment[i])[0], p0[1]);
      stream.lineEnd();
      return;
    }

    subject.push(x = new Intersection(p0, segment, null, true));
    clip.push(x.o = new Intersection(p0, null, x, false));
    subject.push(x = new Intersection(p1, segment, null, false));
    clip.push(x.o = new Intersection(p1, null, x, true));
  });

  if (!subject.length) return;

  clip.sort(compareIntersection);
  polygon_link(subject);
  polygon_link(clip);

  for (i = 0, n = clip.length; i < n; ++i) {
    clip[i].e = startInside = !startInside;
  }

  var start = subject[0],
      points,
      point;

  while (1) {
    // Find first unvisited intersection.
    var current = start,
        isSubject = true;
    while (current.v) if ((current = current.n) === start) return;
    points = current.z;
    stream.lineStart();
    do {
      current.v = current.o.v = true;
      if (current.e) {
        if (isSubject) {
          for (i = 0, n = points.length; i < n; ++i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.n.x, 1, stream);
        }
        current = current.n;
      } else {
        if (isSubject) {
          points = current.p.z;
          for (i = points.length - 1; i >= 0; --i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.p.x, -1, stream);
        }
        current = current.p;
      }
      current = current.o;
      points = current.z;
      isSubject = !isSubject;
    } while (!current.v);
    stream.lineEnd();
  }
});

function polygon_link(array) {
  if (!(n = array.length)) return;
  var n,
      i = 0,
      a = array[0],
      b;
  while (++i < n) {
    a.n = b = array[i];
    b.p = a;
    a = b;
  }
  a.n = b = array[0];
  b.p = a;
}

// EXTERNAL MODULE: ./node_modules/d3-array/src/index.js + 31 modules
var src = __webpack_require__(2);

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/extent.js






var clipMax = 1e9, clipMin = -clipMax;

// TODO Use d3-polygon’s polygonContains here for the ring check?
// TODO Eliminate duplicate buffering in clipBuffer and polygon.push?

function extent_clipExtent(x0, y0, x1, y1) {

  function visible(x, y) {
    return x0 <= x && x <= x1 && y0 <= y && y <= y1;
  }

  function interpolate(from, to, direction, stream) {
    var a = 0, a1 = 0;
    if (from == null
        || (a = corner(from, direction)) !== (a1 = corner(to, direction))
        || comparePoint(from, to) < 0 ^ direction > 0) {
      do stream.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0);
      while ((a = (a + direction + 4) % 4) !== a1);
    } else {
      stream.point(to[0], to[1]);
    }
  }

  function corner(p, direction) {
    return abs(p[0] - x0) < epsilon ? direction > 0 ? 0 : 3
        : abs(p[0] - x1) < epsilon ? direction > 0 ? 2 : 1
        : abs(p[1] - y0) < epsilon ? direction > 0 ? 1 : 0
        : direction > 0 ? 3 : 2; // abs(p[1] - y1) < epsilon
  }

  function compareIntersection(a, b) {
    return comparePoint(a.x, b.x);
  }

  function comparePoint(a, b) {
    var ca = corner(a, 1),
        cb = corner(b, 1);
    return ca !== cb ? ca - cb
        : ca === 0 ? b[1] - a[1]
        : ca === 1 ? a[0] - b[0]
        : ca === 2 ? a[1] - b[1]
        : b[0] - a[0];
  }

  return function(stream) {
    var activeStream = stream,
        bufferStream = buffer(),
        segments,
        polygon,
        ring,
        x__, y__, v__, // first point
        x_, y_, v_, // previous point
        first,
        clean;

    var clipStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: polygonStart,
      polygonEnd: polygonEnd
    };

    function point(x, y) {
      if (visible(x, y)) activeStream.point(x, y);
    }

    function polygonInside() {
      var winding = 0;

      for (var i = 0, n = polygon.length; i < n; ++i) {
        for (var ring = polygon[i], j = 1, m = ring.length, point = ring[0], a0, a1, b0 = point[0], b1 = point[1]; j < m; ++j) {
          a0 = b0, a1 = b1, point = ring[j], b0 = point[0], b1 = point[1];
          if (a1 <= y1) { if (b1 > y1 && (b0 - a0) * (y1 - a1) > (b1 - a1) * (x0 - a0)) ++winding; }
          else { if (b1 <= y1 && (b0 - a0) * (y1 - a1) < (b1 - a1) * (x0 - a0)) --winding; }
        }
      }

      return winding;
    }

    // Buffer geometry within a polygon and then clip it en masse.
    function polygonStart() {
      activeStream = bufferStream, segments = [], polygon = [], clean = true;
    }

    function polygonEnd() {
      var startInside = polygonInside(),
          cleanInside = clean && startInside,
          visible = (segments = Object(src["a" /* merge */])(segments)).length;
      if (cleanInside || visible) {
        stream.polygonStart();
        if (cleanInside) {
          stream.lineStart();
          interpolate(null, null, 1, stream);
          stream.lineEnd();
        }
        if (visible) {
          clip_polygon(segments, compareIntersection, startInside, interpolate, stream);
        }
        stream.polygonEnd();
      }
      activeStream = stream, segments = polygon = ring = null;
    }

    function lineStart() {
      clipStream.point = linePoint;
      if (polygon) polygon.push(ring = []);
      first = true;
      v_ = false;
      x_ = y_ = NaN;
    }

    // TODO rather than special-case polygons, simply handle them separately.
    // Ideally, coincident intersection points should be jittered to avoid
    // clipping issues.
    function lineEnd() {
      if (segments) {
        linePoint(x__, y__);
        if (v__ && v_) bufferStream.rejoin();
        segments.push(bufferStream.result());
      }
      clipStream.point = point;
      if (v_) activeStream.lineEnd();
    }

    function linePoint(x, y) {
      var v = visible(x, y);
      if (polygon) ring.push([x, y]);
      if (first) {
        x__ = x, y__ = y, v__ = v;
        first = false;
        if (v) {
          activeStream.lineStart();
          activeStream.point(x, y);
        }
      } else {
        if (v && v_) activeStream.point(x, y);
        else {
          var a = [x_ = Math.max(clipMin, Math.min(clipMax, x_)), y_ = Math.max(clipMin, Math.min(clipMax, y_))],
              b = [x = Math.max(clipMin, Math.min(clipMax, x)), y = Math.max(clipMin, Math.min(clipMax, y))];
          if (clip_line(a, b, x0, y0, x1, y1)) {
            if (!v_) {
              activeStream.lineStart();
              activeStream.point(a[0], a[1]);
            }
            activeStream.point(b[0], b[1]);
            if (!v) activeStream.lineEnd();
            clean = false;
          } else if (v) {
            activeStream.lineStart();
            activeStream.point(x, y);
            clean = false;
          }
        }
      }
      x_ = x, y_ = y, v_ = v;
    }

    return clipStream;
  };
}

/* harmony default export */ var clip_extent = (function() {
  var x0 = 0,
      y0 = 0,
      x1 = 960,
      y1 = 500,
      cache,
      cacheStream,
      clip;

  return clip = {
    stream: function(stream) {
      return cache && cacheStream === stream ? cache : cache = extent_clipExtent(x0, y0, x1, y1)(cacheStream = stream);
    },
    extent: function(_) {
      return arguments.length ? (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1], cache = cacheStream = null, clip) : [[x0, y0], [x1, y1]];
    }
  };
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/polygonContains.js




var sum = adder();

/* harmony default export */ var polygonContains = (function(polygon, point) {
  var lambda = point[0],
      phi = point[1],
      normal = [sin(lambda), -cos(lambda), 0],
      angle = 0,
      winding = 0;

  sum.reset();

  for (var i = 0, n = polygon.length; i < n; ++i) {
    if (!(m = (ring = polygon[i]).length)) continue;
    var ring,
        m,
        point0 = ring[m - 1],
        lambda0 = point0[0],
        phi0 = point0[1] / 2 + quarterPi,
        sinPhi0 = sin(phi0),
        cosPhi0 = cos(phi0);

    for (var j = 0; j < m; ++j, lambda0 = lambda1, sinPhi0 = sinPhi1, cosPhi0 = cosPhi1, point0 = point1) {
      var point1 = ring[j],
          lambda1 = point1[0],
          phi1 = point1[1] / 2 + quarterPi,
          sinPhi1 = sin(phi1),
          cosPhi1 = cos(phi1),
          delta = lambda1 - lambda0,
          sign = delta >= 0 ? 1 : -1,
          absDelta = sign * delta,
          antimeridian = absDelta > pi,
          k = sinPhi0 * sinPhi1;

      sum.add(atan2(k * sign * sin(absDelta), cosPhi0 * cosPhi1 + k * cos(absDelta)));
      angle += antimeridian ? delta + sign * tau : delta;

      // Are the longitudes either side of the point’s meridian (lambda),
      // and are the latitudes smaller than the parallel (phi)?
      if (antimeridian ^ lambda0 >= lambda ^ lambda1 >= lambda) {
        var arc = cartesianCross(cartesian_cartesian(point0), cartesian_cartesian(point1));
        cartesianNormalizeInPlace(arc);
        var intersection = cartesianCross(normal, arc);
        cartesianNormalizeInPlace(intersection);
        var phiArc = (antimeridian ^ delta >= 0 ? -1 : 1) * asin(intersection[2]);
        if (phi > phiArc || phi === phiArc && (arc[0] || arc[1])) {
          winding += antimeridian ^ delta >= 0 ? 1 : -1;
        }
      }
    }
  }

  // First, determine whether the South pole is inside or outside:
  //
  // It is inside if:
  // * the polygon winds around it in a clockwise direction.
  // * the polygon does not (cumulatively) wind around it, but has a negative
  //   (counter-clockwise) area.
  //
  // Second, count the (signed) number of times a segment crosses a lambda
  // from the point to the South pole.  If it is zero, then the point is the
  // same side as the South pole.

  return (angle < -epsilon || angle < epsilon && sum < -epsilon) ^ (winding & 1);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/length.js





var lengthSum = adder(),
    length_lambda0,
    length_sinPhi0,
    length_cosPhi0;

var lengthStream = {
  sphere: noop,
  point: noop,
  lineStart: lengthLineStart,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop
};

function lengthLineStart() {
  lengthStream.point = lengthPointFirst;
  lengthStream.lineEnd = lengthLineEnd;
}

function lengthLineEnd() {
  lengthStream.point = lengthStream.lineEnd = noop;
}

function lengthPointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  length_lambda0 = lambda, length_sinPhi0 = sin(phi), length_cosPhi0 = cos(phi);
  lengthStream.point = lengthPoint;
}

function lengthPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var sinPhi = sin(phi),
      cosPhi = cos(phi),
      delta = abs(lambda - length_lambda0),
      cosDelta = cos(delta),
      sinDelta = sin(delta),
      x = cosPhi * sinDelta,
      y = length_cosPhi0 * sinPhi - length_sinPhi0 * cosPhi * cosDelta,
      z = length_sinPhi0 * sinPhi + length_cosPhi0 * cosPhi * cosDelta;
  lengthSum.add(atan2(sqrt(x * x + y * y), z));
  length_lambda0 = lambda, length_sinPhi0 = sinPhi, length_cosPhi0 = cosPhi;
}

/* harmony default export */ var src_length = (function(object) {
  lengthSum.reset();
  src_stream(object, lengthStream);
  return +lengthSum;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/distance.js


var distance_coordinates = [null, null],
    distance_object = {type: "LineString", coordinates: distance_coordinates};

/* harmony default export */ var distance = (function(a, b) {
  distance_coordinates[0] = a;
  distance_coordinates[1] = b;
  return src_length(distance_object);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/contains.js




var containsObjectType = {
  Feature: function(object, point) {
    return containsGeometry(object.geometry, point);
  },
  FeatureCollection: function(object, point) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n) if (containsGeometry(features[i].geometry, point)) return true;
    return false;
  }
};

var containsGeometryType = {
  Sphere: function() {
    return true;
  },
  Point: function(object, point) {
    return containsPoint(object.coordinates, point);
  },
  MultiPoint: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsPoint(coordinates[i], point)) return true;
    return false;
  },
  LineString: function(object, point) {
    return containsLine(object.coordinates, point);
  },
  MultiLineString: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsLine(coordinates[i], point)) return true;
    return false;
  },
  Polygon: function(object, point) {
    return containsPolygon(object.coordinates, point);
  },
  MultiPolygon: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsPolygon(coordinates[i], point)) return true;
    return false;
  },
  GeometryCollection: function(object, point) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n) if (containsGeometry(geometries[i], point)) return true;
    return false;
  }
};

function containsGeometry(geometry, point) {
  return geometry && containsGeometryType.hasOwnProperty(geometry.type)
      ? containsGeometryType[geometry.type](geometry, point)
      : false;
}

function containsPoint(coordinates, point) {
  return distance(coordinates, point) === 0;
}

function containsLine(coordinates, point) {
  var ab = distance(coordinates[0], coordinates[1]),
      ao = distance(coordinates[0], point),
      ob = distance(point, coordinates[1]);
  return ao + ob <= ab + epsilon;
}

function containsPolygon(coordinates, point) {
  return !!polygonContains(coordinates.map(ringRadians), pointRadians(point));
}

function ringRadians(ring) {
  return ring = ring.map(pointRadians), ring.pop(), ring;
}

function pointRadians(point) {
  return [point[0] * radians, point[1] * radians];
}

/* harmony default export */ var contains = (function(object, point) {
  return (object && containsObjectType.hasOwnProperty(object.type)
      ? containsObjectType[object.type]
      : containsGeometry)(object, point);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/graticule.js



function graticuleX(y0, y1, dy) {
  var y = Object(src["b" /* range */])(y0, y1 - epsilon, dy).concat(y1);
  return function(x) { return y.map(function(y) { return [x, y]; }); };
}

function graticuleY(x0, x1, dx) {
  var x = Object(src["b" /* range */])(x0, x1 - epsilon, dx).concat(x1);
  return function(y) { return x.map(function(x) { return [x, y]; }); };
}

function graticule_graticule() {
  var x1, x0, X1, X0,
      y1, y0, Y1, Y0,
      dx = 10, dy = dx, DX = 90, DY = 360,
      x, y, X, Y,
      precision = 2.5;

  function graticule() {
    return {type: "MultiLineString", coordinates: lines()};
  }

  function lines() {
    return Object(src["b" /* range */])(ceil(X0 / DX) * DX, X1, DX).map(X)
        .concat(Object(src["b" /* range */])(ceil(Y0 / DY) * DY, Y1, DY).map(Y))
        .concat(Object(src["b" /* range */])(ceil(x0 / dx) * dx, x1, dx).filter(function(x) { return abs(x % DX) > epsilon; }).map(x))
        .concat(Object(src["b" /* range */])(ceil(y0 / dy) * dy, y1, dy).filter(function(y) { return abs(y % DY) > epsilon; }).map(y));
  }

  graticule.lines = function() {
    return lines().map(function(coordinates) { return {type: "LineString", coordinates: coordinates}; });
  };

  graticule.outline = function() {
    return {
      type: "Polygon",
      coordinates: [
        X(X0).concat(
        Y(Y1).slice(1),
        X(X1).reverse().slice(1),
        Y(Y0).reverse().slice(1))
      ]
    };
  };

  graticule.extent = function(_) {
    if (!arguments.length) return graticule.extentMinor();
    return graticule.extentMajor(_).extentMinor(_);
  };

  graticule.extentMajor = function(_) {
    if (!arguments.length) return [[X0, Y0], [X1, Y1]];
    X0 = +_[0][0], X1 = +_[1][0];
    Y0 = +_[0][1], Y1 = +_[1][1];
    if (X0 > X1) _ = X0, X0 = X1, X1 = _;
    if (Y0 > Y1) _ = Y0, Y0 = Y1, Y1 = _;
    return graticule.precision(precision);
  };

  graticule.extentMinor = function(_) {
    if (!arguments.length) return [[x0, y0], [x1, y1]];
    x0 = +_[0][0], x1 = +_[1][0];
    y0 = +_[0][1], y1 = +_[1][1];
    if (x0 > x1) _ = x0, x0 = x1, x1 = _;
    if (y0 > y1) _ = y0, y0 = y1, y1 = _;
    return graticule.precision(precision);
  };

  graticule.step = function(_) {
    if (!arguments.length) return graticule.stepMinor();
    return graticule.stepMajor(_).stepMinor(_);
  };

  graticule.stepMajor = function(_) {
    if (!arguments.length) return [DX, DY];
    DX = +_[0], DY = +_[1];
    return graticule;
  };

  graticule.stepMinor = function(_) {
    if (!arguments.length) return [dx, dy];
    dx = +_[0], dy = +_[1];
    return graticule;
  };

  graticule.precision = function(_) {
    if (!arguments.length) return precision;
    precision = +_;
    x = graticuleX(y0, y1, 90);
    y = graticuleY(x0, x1, precision);
    X = graticuleX(Y0, Y1, 90);
    Y = graticuleY(X0, X1, precision);
    return graticule;
  };

  return graticule
      .extentMajor([[-180, -90 + epsilon], [180, 90 - epsilon]])
      .extentMinor([[-180, -80 - epsilon], [180, 80 + epsilon]]);
}

function graticule10() {
  return graticule_graticule()();
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/interpolate.js


/* harmony default export */ var src_interpolate = (function(a, b) {
  var x0 = a[0] * radians,
      y0 = a[1] * radians,
      x1 = b[0] * radians,
      y1 = b[1] * radians,
      cy0 = cos(y0),
      sy0 = sin(y0),
      cy1 = cos(y1),
      sy1 = sin(y1),
      kx0 = cy0 * cos(x0),
      ky0 = cy0 * sin(x0),
      kx1 = cy1 * cos(x1),
      ky1 = cy1 * sin(x1),
      d = 2 * asin(sqrt(haversin(y1 - y0) + cy0 * cy1 * haversin(x1 - x0))),
      k = sin(d);

  var interpolate = d ? function(t) {
    var B = sin(t *= d) / k,
        A = sin(d - t) / k,
        x = A * kx0 + B * kx1,
        y = A * ky0 + B * ky1,
        z = A * sy0 + B * sy1;
    return [
      atan2(y, x) * degrees,
      atan2(z, sqrt(x * x + y * y)) * degrees
    ];
  } : function() {
    return [x0 * degrees, y0 * degrees];
  };

  interpolate.distance = d;

  return interpolate;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/identity.js
/* harmony default export */ var identity = (function(x) {
  return x;
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/area.js




var area_areaSum = adder(),
    area_areaRingSum = adder(),
    area_x00,
    area_y00,
    area_x0,
    area_y0;

var area_areaStream = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    area_areaStream.lineStart = area_areaRingStart;
    area_areaStream.lineEnd = area_areaRingEnd;
  },
  polygonEnd: function() {
    area_areaStream.lineStart = area_areaStream.lineEnd = area_areaStream.point = noop;
    area_areaSum.add(abs(area_areaRingSum));
    area_areaRingSum.reset();
  },
  result: function() {
    var area = area_areaSum / 2;
    area_areaSum.reset();
    return area;
  }
};

function area_areaRingStart() {
  area_areaStream.point = area_areaPointFirst;
}

function area_areaPointFirst(x, y) {
  area_areaStream.point = area_areaPoint;
  area_x00 = area_x0 = x, area_y00 = area_y0 = y;
}

function area_areaPoint(x, y) {
  area_areaRingSum.add(area_y0 * x - area_x0 * y);
  area_x0 = x, area_y0 = y;
}

function area_areaRingEnd() {
  area_areaPoint(area_x00, area_y00);
}

/* harmony default export */ var path_area = (area_areaStream);

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/bounds.js


var bounds_x0 = Infinity,
    bounds_y0 = bounds_x0,
    bounds_x1 = -bounds_x0,
    bounds_y1 = bounds_x1;

var bounds_boundsStream = {
  point: bounds_boundsPoint,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop,
  result: function() {
    var bounds = [[bounds_x0, bounds_y0], [bounds_x1, bounds_y1]];
    bounds_x1 = bounds_y1 = -(bounds_y0 = bounds_x0 = Infinity);
    return bounds;
  }
};

function bounds_boundsPoint(x, y) {
  if (x < bounds_x0) bounds_x0 = x;
  if (x > bounds_x1) bounds_x1 = x;
  if (y < bounds_y0) bounds_y0 = y;
  if (y > bounds_y1) bounds_y1 = y;
}

/* harmony default export */ var path_bounds = (bounds_boundsStream);

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/centroid.js


// TODO Enforce positive area for exterior, negative area for interior?

var path_centroid_X0 = 0,
    path_centroid_Y0 = 0,
    centroid_Z0 = 0,
    path_centroid_X1 = 0,
    path_centroid_Y1 = 0,
    centroid_Z1 = 0,
    centroid_X2 = 0,
    centroid_Y2 = 0,
    centroid_Z2 = 0,
    centroid_x00,
    centroid_y00,
    path_centroid_x0,
    path_centroid_y0;

var centroid_centroidStream = {
  point: centroid_centroidPoint,
  lineStart: centroid_centroidLineStart,
  lineEnd: centroid_centroidLineEnd,
  polygonStart: function() {
    centroid_centroidStream.lineStart = centroid_centroidRingStart;
    centroid_centroidStream.lineEnd = centroid_centroidRingEnd;
  },
  polygonEnd: function() {
    centroid_centroidStream.point = centroid_centroidPoint;
    centroid_centroidStream.lineStart = centroid_centroidLineStart;
    centroid_centroidStream.lineEnd = centroid_centroidLineEnd;
  },
  result: function() {
    var centroid = centroid_Z2 ? [centroid_X2 / centroid_Z2, centroid_Y2 / centroid_Z2]
        : centroid_Z1 ? [path_centroid_X1 / centroid_Z1, path_centroid_Y1 / centroid_Z1]
        : centroid_Z0 ? [path_centroid_X0 / centroid_Z0, path_centroid_Y0 / centroid_Z0]
        : [NaN, NaN];
    path_centroid_X0 = path_centroid_Y0 = centroid_Z0 =
    path_centroid_X1 = path_centroid_Y1 = centroid_Z1 =
    centroid_X2 = centroid_Y2 = centroid_Z2 = 0;
    return centroid;
  }
};

function centroid_centroidPoint(x, y) {
  path_centroid_X0 += x;
  path_centroid_Y0 += y;
  ++centroid_Z0;
}

function centroid_centroidLineStart() {
  centroid_centroidStream.point = centroidPointFirstLine;
}

function centroidPointFirstLine(x, y) {
  centroid_centroidStream.point = centroidPointLine;
  centroid_centroidPoint(path_centroid_x0 = x, path_centroid_y0 = y);
}

function centroidPointLine(x, y) {
  var dx = x - path_centroid_x0, dy = y - path_centroid_y0, z = sqrt(dx * dx + dy * dy);
  path_centroid_X1 += z * (path_centroid_x0 + x) / 2;
  path_centroid_Y1 += z * (path_centroid_y0 + y) / 2;
  centroid_Z1 += z;
  centroid_centroidPoint(path_centroid_x0 = x, path_centroid_y0 = y);
}

function centroid_centroidLineEnd() {
  centroid_centroidStream.point = centroid_centroidPoint;
}

function centroid_centroidRingStart() {
  centroid_centroidStream.point = centroidPointFirstRing;
}

function centroid_centroidRingEnd() {
  centroidPointRing(centroid_x00, centroid_y00);
}

function centroidPointFirstRing(x, y) {
  centroid_centroidStream.point = centroidPointRing;
  centroid_centroidPoint(centroid_x00 = path_centroid_x0 = x, centroid_y00 = path_centroid_y0 = y);
}

function centroidPointRing(x, y) {
  var dx = x - path_centroid_x0,
      dy = y - path_centroid_y0,
      z = sqrt(dx * dx + dy * dy);

  path_centroid_X1 += z * (path_centroid_x0 + x) / 2;
  path_centroid_Y1 += z * (path_centroid_y0 + y) / 2;
  centroid_Z1 += z;

  z = path_centroid_y0 * x - path_centroid_x0 * y;
  centroid_X2 += z * (path_centroid_x0 + x);
  centroid_Y2 += z * (path_centroid_y0 + y);
  centroid_Z2 += z * 3;
  centroid_centroidPoint(path_centroid_x0 = x, path_centroid_y0 = y);
}

/* harmony default export */ var path_centroid = (centroid_centroidStream);

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/context.js



function PathContext(context) {
  this._context = context;
}

PathContext.prototype = {
  _radius: 4.5,
  pointRadius: function(_) {
    return this._radius = _, this;
  },
  polygonStart: function() {
    this._line = 0;
  },
  polygonEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line === 0) this._context.closePath();
    this._point = NaN;
  },
  point: function(x, y) {
    switch (this._point) {
      case 0: {
        this._context.moveTo(x, y);
        this._point = 1;
        break;
      }
      case 1: {
        this._context.lineTo(x, y);
        break;
      }
      default: {
        this._context.moveTo(x + this._radius, y);
        this._context.arc(x, y, this._radius, 0, tau);
        break;
      }
    }
  },
  result: noop
};

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/measure.js




var measure_lengthSum = adder(),
    lengthRing,
    measure_x00,
    measure_y00,
    measure_x0,
    measure_y0;

var measure_lengthStream = {
  point: noop,
  lineStart: function() {
    measure_lengthStream.point = measure_lengthPointFirst;
  },
  lineEnd: function() {
    if (lengthRing) measure_lengthPoint(measure_x00, measure_y00);
    measure_lengthStream.point = noop;
  },
  polygonStart: function() {
    lengthRing = true;
  },
  polygonEnd: function() {
    lengthRing = null;
  },
  result: function() {
    var length = +measure_lengthSum;
    measure_lengthSum.reset();
    return length;
  }
};

function measure_lengthPointFirst(x, y) {
  measure_lengthStream.point = measure_lengthPoint;
  measure_x00 = measure_x0 = x, measure_y00 = measure_y0 = y;
}

function measure_lengthPoint(x, y) {
  measure_x0 -= x, measure_y0 -= y;
  measure_lengthSum.add(sqrt(measure_x0 * measure_x0 + measure_y0 * measure_y0));
  measure_x0 = x, measure_y0 = y;
}

/* harmony default export */ var measure = (measure_lengthStream);

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/string.js
function PathString() {
  this._string = [];
}

PathString.prototype = {
  _circle: string_circle(4.5),
  pointRadius: function(_) {
    return this._circle = string_circle(_), this;
  },
  polygonStart: function() {
    this._line = 0;
  },
  polygonEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line === 0) this._string.push("Z");
    this._point = NaN;
  },
  point: function(x, y) {
    switch (this._point) {
      case 0: {
        this._string.push("M", x, ",", y);
        this._point = 1;
        break;
      }
      case 1: {
        this._string.push("L", x, ",", y);
        break;
      }
      default: {
        this._string.push("M", x, ",", y, this._circle);
        break;
      }
    }
  },
  result: function() {
    if (this._string.length) {
      var result = this._string.join("");
      this._string = [];
      return result;
    }
  }
};

function string_circle(radius) {
  return "m0," + radius
      + "a" + radius + "," + radius + " 0 1,1 0," + -2 * radius
      + "a" + radius + "," + radius + " 0 1,1 0," + 2 * radius
      + "z";
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/path/index.js









/* harmony default export */ var src_path = (function(projection, context) {
  var pointRadius = 4.5,
      projectionStream,
      contextStream;

  function path(object) {
    if (object) {
      if (typeof pointRadius === "function") contextStream.pointRadius(+pointRadius.apply(this, arguments));
      src_stream(object, projectionStream(contextStream));
    }
    return contextStream.result();
  }

  path.area = function(object) {
    src_stream(object, projectionStream(path_area));
    return path_area.result();
  };

  path.measure = function(object) {
    src_stream(object, projectionStream(measure));
    return measure.result();
  };

  path.bounds = function(object) {
    src_stream(object, projectionStream(path_bounds));
    return path_bounds.result();
  };

  path.centroid = function(object) {
    src_stream(object, projectionStream(path_centroid));
    return path_centroid.result();
  };

  path.projection = function(_) {
    return arguments.length ? (projectionStream = _ == null ? (projection = null, identity) : (projection = _).stream, path) : projection;
  };

  path.context = function(_) {
    if (!arguments.length) return context;
    contextStream = _ == null ? (context = null, new PathString) : new PathContext(context = _);
    if (typeof pointRadius !== "function") contextStream.pointRadius(pointRadius);
    return path;
  };

  path.pointRadius = function(_) {
    if (!arguments.length) return pointRadius;
    pointRadius = typeof _ === "function" ? _ : (contextStream.pointRadius(+_), +_);
    return path;
  };

  return path.projection(projection).context(context);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/index.js






/* harmony default export */ var src_clip = (function(pointVisible, clipLine, interpolate, start) {
  return function(rotate, sink) {
    var line = clipLine(sink),
        rotatedStart = rotate.invert(start[0], start[1]),
        ringBuffer = buffer(),
        ringSink = clipLine(ringBuffer),
        polygonStarted = false,
        polygon,
        segments,
        ring;

    var clip = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function() {
        clip.point = pointRing;
        clip.lineStart = ringStart;
        clip.lineEnd = ringEnd;
        segments = [];
        polygon = [];
      },
      polygonEnd: function() {
        clip.point = point;
        clip.lineStart = lineStart;
        clip.lineEnd = lineEnd;
        segments = Object(src["a" /* merge */])(segments);
        var startInside = polygonContains(polygon, rotatedStart);
        if (segments.length) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          clip_polygon(segments, clip_compareIntersection, startInside, interpolate, sink);
        } else if (startInside) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          interpolate(null, null, 1, sink);
          sink.lineEnd();
        }
        if (polygonStarted) sink.polygonEnd(), polygonStarted = false;
        segments = polygon = null;
      },
      sphere: function() {
        sink.polygonStart();
        sink.lineStart();
        interpolate(null, null, 1, sink);
        sink.lineEnd();
        sink.polygonEnd();
      }
    };

    function point(lambda, phi) {
      var point = rotate(lambda, phi);
      if (pointVisible(lambda = point[0], phi = point[1])) sink.point(lambda, phi);
    }

    function pointLine(lambda, phi) {
      var point = rotate(lambda, phi);
      line.point(point[0], point[1]);
    }

    function lineStart() {
      clip.point = pointLine;
      line.lineStart();
    }

    function lineEnd() {
      clip.point = point;
      line.lineEnd();
    }

    function pointRing(lambda, phi) {
      ring.push([lambda, phi]);
      var point = rotate(lambda, phi);
      ringSink.point(point[0], point[1]);
    }

    function ringStart() {
      ringSink.lineStart();
      ring = [];
    }

    function ringEnd() {
      pointRing(ring[0][0], ring[0][1]);
      ringSink.lineEnd();

      var clean = ringSink.clean(),
          ringSegments = ringBuffer.result(),
          i, n = ringSegments.length, m,
          segment,
          point;

      ring.pop();
      polygon.push(ring);
      ring = null;

      if (!n) return;

      // No intersections.
      if (clean & 1) {
        segment = ringSegments[0];
        if ((m = segment.length - 1) > 0) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          for (i = 0; i < m; ++i) sink.point((point = segment[i])[0], point[1]);
          sink.lineEnd();
        }
        return;
      }

      // Rejoin connected segments.
      // TODO reuse ringBuffer.rejoin()?
      if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));

      segments.push(ringSegments.filter(validSegment));
    }

    return clip;
  };
});

function validSegment(segment) {
  return segment.length > 1;
}

// Intersections are sorted along the clip edge. For both antimeridian cutting
// and circle clipping, the same comparison is used.
function clip_compareIntersection(a, b) {
  return ((a = a.x)[0] < 0 ? a[1] - halfPi - epsilon : halfPi - a[1])
       - ((b = b.x)[0] < 0 ? b[1] - halfPi - epsilon : halfPi - b[1]);
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/antimeridian.js



/* harmony default export */ var clip_antimeridian = (src_clip(
  function() { return true; },
  clipAntimeridianLine,
  clipAntimeridianInterpolate,
  [-pi, -halfPi]
));

// Takes a line and cuts into visible segments. Return values: 0 - there were
// intersections or the line was empty; 1 - no intersections; 2 - there were
// intersections, and the first and last segments should be rejoined.
function clipAntimeridianLine(stream) {
  var lambda0 = NaN,
      phi0 = NaN,
      sign0 = NaN,
      clean; // no intersections

  return {
    lineStart: function() {
      stream.lineStart();
      clean = 1;
    },
    point: function(lambda1, phi1) {
      var sign1 = lambda1 > 0 ? pi : -pi,
          delta = abs(lambda1 - lambda0);
      if (abs(delta - pi) < epsilon) { // line crosses a pole
        stream.point(lambda0, phi0 = (phi0 + phi1) / 2 > 0 ? halfPi : -halfPi);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        stream.point(lambda1, phi0);
        clean = 0;
      } else if (sign0 !== sign1 && delta >= pi) { // line crosses antimeridian
        if (abs(lambda0 - sign0) < epsilon) lambda0 -= sign0 * epsilon; // handle degeneracies
        if (abs(lambda1 - sign1) < epsilon) lambda1 -= sign1 * epsilon;
        phi0 = clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        clean = 0;
      }
      stream.point(lambda0 = lambda1, phi0 = phi1);
      sign0 = sign1;
    },
    lineEnd: function() {
      stream.lineEnd();
      lambda0 = phi0 = NaN;
    },
    clean: function() {
      return 2 - clean; // if intersections, rejoin first and last segments
    }
  };
}

function clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1) {
  var cosPhi0,
      cosPhi1,
      sinLambda0Lambda1 = sin(lambda0 - lambda1);
  return abs(sinLambda0Lambda1) > epsilon
      ? atan((sin(phi0) * (cosPhi1 = cos(phi1)) * sin(lambda1)
          - sin(phi1) * (cosPhi0 = cos(phi0)) * sin(lambda0))
          / (cosPhi0 * cosPhi1 * sinLambda0Lambda1))
      : (phi0 + phi1) / 2;
}

function clipAntimeridianInterpolate(from, to, direction, stream) {
  var phi;
  if (from == null) {
    phi = direction * halfPi;
    stream.point(-pi, phi);
    stream.point(0, phi);
    stream.point(pi, phi);
    stream.point(pi, 0);
    stream.point(pi, -phi);
    stream.point(0, -phi);
    stream.point(-pi, -phi);
    stream.point(-pi, 0);
    stream.point(-pi, phi);
  } else if (abs(from[0] - to[0]) > epsilon) {
    var lambda = from[0] < to[0] ? pi : -pi;
    phi = direction * lambda / 2;
    stream.point(-lambda, phi);
    stream.point(0, phi);
    stream.point(lambda, phi);
  } else {
    stream.point(to[0], to[1]);
  }
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/clip/circle.js






/* harmony default export */ var clip_circle = (function(radius, delta) {
  var cr = cos(radius),
      smallRadius = cr > 0,
      notHemisphere = abs(cr) > epsilon; // TODO optimise for this common case

  function interpolate(from, to, direction, stream) {
    circleStream(stream, radius, delta, direction, from, to);
  }

  function visible(lambda, phi) {
    return cos(lambda) * cos(phi) > cr;
  }

  // Takes a line and cuts into visible segments. Return values used for polygon
  // clipping: 0 - there were intersections or the line was empty; 1 - no
  // intersections 2 - there were intersections, and the first and last segments
  // should be rejoined.
  function clipLine(stream) {
    var point0, // previous point
        c0, // code for previous point
        v0, // visibility of previous point
        v00, // visibility of first point
        clean; // no intersections
    return {
      lineStart: function() {
        v00 = v0 = false;
        clean = 1;
      },
      point: function(lambda, phi) {
        var point1 = [lambda, phi],
            point2,
            v = visible(lambda, phi),
            c = smallRadius
              ? v ? 0 : code(lambda, phi)
              : v ? code(lambda + (lambda < 0 ? pi : -pi), phi) : 0;
        if (!point0 && (v00 = v0 = v)) stream.lineStart();
        // Handle degeneracies.
        // TODO ignore if not clipping polygons.
        if (v !== v0) {
          point2 = intersect(point0, point1);
          if (pointEqual(point0, point2) || pointEqual(point1, point2)) {
            point1[0] += epsilon;
            point1[1] += epsilon;
            v = visible(point1[0], point1[1]);
          }
        }
        if (v !== v0) {
          clean = 0;
          if (v) {
            // outside going in
            stream.lineStart();
            point2 = intersect(point1, point0);
            stream.point(point2[0], point2[1]);
          } else {
            // inside going out
            point2 = intersect(point0, point1);
            stream.point(point2[0], point2[1]);
            stream.lineEnd();
          }
          point0 = point2;
        } else if (notHemisphere && point0 && smallRadius ^ v) {
          var t;
          // If the codes for two points are different, or are both zero,
          // and there this segment intersects with the small circle.
          if (!(c & c0) && (t = intersect(point1, point0, true))) {
            clean = 0;
            if (smallRadius) {
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
            } else {
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
            }
          }
        }
        if (v && (!point0 || !pointEqual(point0, point1))) {
          stream.point(point1[0], point1[1]);
        }
        point0 = point1, v0 = v, c0 = c;
      },
      lineEnd: function() {
        if (v0) stream.lineEnd();
        point0 = null;
      },
      // Rejoin first and last segments if there were intersections and the first
      // and last points were visible.
      clean: function() {
        return clean | ((v00 && v0) << 1);
      }
    };
  }

  // Intersects the great circle between a and b with the clip circle.
  function intersect(a, b, two) {
    var pa = cartesian_cartesian(a),
        pb = cartesian_cartesian(b);

    // We have two planes, n1.p = d1 and n2.p = d2.
    // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).
    var n1 = [1, 0, 0], // normal
        n2 = cartesianCross(pa, pb),
        n2n2 = cartesianDot(n2, n2),
        n1n2 = n2[0], // cartesianDot(n1, n2),
        determinant = n2n2 - n1n2 * n1n2;

    // Two polar points.
    if (!determinant) return !two && a;

    var c1 =  cr * n2n2 / determinant,
        c2 = -cr * n1n2 / determinant,
        n1xn2 = cartesianCross(n1, n2),
        A = cartesianScale(n1, c1),
        B = cartesianScale(n2, c2);
    cartesianAddInPlace(A, B);

    // Solve |p(t)|^2 = 1.
    var u = n1xn2,
        w = cartesianDot(A, u),
        uu = cartesianDot(u, u),
        t2 = w * w - uu * (cartesianDot(A, A) - 1);

    if (t2 < 0) return;

    var t = sqrt(t2),
        q = cartesianScale(u, (-w - t) / uu);
    cartesianAddInPlace(q, A);
    q = cartesian_spherical(q);

    if (!two) return q;

    // Two intersection points.
    var lambda0 = a[0],
        lambda1 = b[0],
        phi0 = a[1],
        phi1 = b[1],
        z;

    if (lambda1 < lambda0) z = lambda0, lambda0 = lambda1, lambda1 = z;

    var delta = lambda1 - lambda0,
        polar = abs(delta - pi) < epsilon,
        meridian = polar || delta < epsilon;

    if (!polar && phi1 < phi0) z = phi0, phi0 = phi1, phi1 = z;

    // Check that the first point is between a and b.
    if (meridian
        ? polar
          ? phi0 + phi1 > 0 ^ q[1] < (abs(q[0] - lambda0) < epsilon ? phi0 : phi1)
          : phi0 <= q[1] && q[1] <= phi1
        : delta > pi ^ (lambda0 <= q[0] && q[0] <= lambda1)) {
      var q1 = cartesianScale(u, (-w + t) / uu);
      cartesianAddInPlace(q1, A);
      return [q, cartesian_spherical(q1)];
    }
  }

  // Generates a 4-bit vector representing the location of a point relative to
  // the small circle's bounding box.
  function code(lambda, phi) {
    var r = smallRadius ? radius : pi - radius,
        code = 0;
    if (lambda < -r) code |= 1; // left
    else if (lambda > r) code |= 2; // right
    if (phi < -r) code |= 4; // below
    else if (phi > r) code |= 8; // above
    return code;
  }

  return src_clip(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-pi, radius - pi]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/transform.js
/* harmony default export */ var src_transform = (function(methods) {
  return {
    stream: transformer(methods)
  };
});

function transformer(methods) {
  return function(stream) {
    var s = new TransformStream;
    for (var key in methods) s[key] = methods[key];
    s.stream = stream;
    return s;
  };
}

function TransformStream() {}

TransformStream.prototype = {
  constructor: TransformStream,
  point: function(x, y) { this.stream.point(x, y); },
  sphere: function() { this.stream.sphere(); },
  lineStart: function() { this.stream.lineStart(); },
  lineEnd: function() { this.stream.lineEnd(); },
  polygonStart: function() { this.stream.polygonStart(); },
  polygonEnd: function() { this.stream.polygonEnd(); }
};

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/fit.js



function fitExtent(projection, extent, object) {
  var w = extent[1][0] - extent[0][0],
      h = extent[1][1] - extent[0][1],
      clip = projection.clipExtent && projection.clipExtent();

  projection
      .scale(150)
      .translate([0, 0]);

  if (clip != null) projection.clipExtent(null);

  src_stream(object, projection.stream(path_bounds));

  var b = path_bounds.result(),
      k = Math.min(w / (b[1][0] - b[0][0]), h / (b[1][1] - b[0][1])),
      x = +extent[0][0] + (w - k * (b[1][0] + b[0][0])) / 2,
      y = +extent[0][1] + (h - k * (b[1][1] + b[0][1])) / 2;

  if (clip != null) projection.clipExtent(clip);

  return projection
      .scale(k * 150)
      .translate([x, y]);
}

function fitSize(projection, size, object) {
  return fitExtent(projection, [[0, 0], size], object);
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/resample.js




var maxDepth = 16, // maximum depth of subdivision
    cosMinDistance = cos(30 * radians); // cos(minimum angular distance)

/* harmony default export */ var resample = (function(project, delta2) {
  return +delta2 ? resample_resample(project, delta2) : resampleNone(project);
});

function resampleNone(project) {
  return transformer({
    point: function(x, y) {
      x = project(x, y);
      this.stream.point(x[0], x[1]);
    }
  });
}

function resample_resample(project, delta2) {

  function resampleLineTo(x0, y0, lambda0, a0, b0, c0, x1, y1, lambda1, a1, b1, c1, depth, stream) {
    var dx = x1 - x0,
        dy = y1 - y0,
        d2 = dx * dx + dy * dy;
    if (d2 > 4 * delta2 && depth--) {
      var a = a0 + a1,
          b = b0 + b1,
          c = c0 + c1,
          m = sqrt(a * a + b * b + c * c),
          phi2 = asin(c /= m),
          lambda2 = abs(abs(c) - 1) < epsilon || abs(lambda0 - lambda1) < epsilon ? (lambda0 + lambda1) / 2 : atan2(b, a),
          p = project(lambda2, phi2),
          x2 = p[0],
          y2 = p[1],
          dx2 = x2 - x0,
          dy2 = y2 - y0,
          dz = dy * dx2 - dx * dy2;
      if (dz * dz / d2 > delta2 // perpendicular projected distance
          || abs((dx * dx2 + dy * dy2) / d2 - 0.5) > 0.3 // midpoint close to an end
          || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) { // angular distance
        resampleLineTo(x0, y0, lambda0, a0, b0, c0, x2, y2, lambda2, a /= m, b /= m, c, depth, stream);
        stream.point(x2, y2);
        resampleLineTo(x2, y2, lambda2, a, b, c, x1, y1, lambda1, a1, b1, c1, depth, stream);
      }
    }
  }
  return function(stream) {
    var lambda00, x00, y00, a00, b00, c00, // first point
        lambda0, x0, y0, a0, b0, c0; // previous point

    var resampleStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function() { stream.polygonStart(); resampleStream.lineStart = ringStart; },
      polygonEnd: function() { stream.polygonEnd(); resampleStream.lineStart = lineStart; }
    };

    function point(x, y) {
      x = project(x, y);
      stream.point(x[0], x[1]);
    }

    function lineStart() {
      x0 = NaN;
      resampleStream.point = linePoint;
      stream.lineStart();
    }

    function linePoint(lambda, phi) {
      var c = cartesian_cartesian([lambda, phi]), p = project(lambda, phi);
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x0 = p[0], y0 = p[1], lambda0 = lambda, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
      stream.point(x0, y0);
    }

    function lineEnd() {
      resampleStream.point = point;
      stream.lineEnd();
    }

    function ringStart() {
      lineStart();
      resampleStream.point = ringPoint;
      resampleStream.lineEnd = ringEnd;
    }

    function ringPoint(lambda, phi) {
      linePoint(lambda00 = lambda, phi), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;
      resampleStream.point = linePoint;
    }

    function ringEnd() {
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x00, y00, lambda00, a00, b00, c00, maxDepth, stream);
      resampleStream.lineEnd = lineEnd;
      lineEnd();
    }

    return resampleStream;
  };
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/index.js











var transformRadians = transformer({
  point: function(x, y) {
    this.stream.point(x * radians, y * radians);
  }
});

function projection_projection(project) {
  return projectionMutator(function() { return project; })();
}

function projectionMutator(projectAt) {
  var project,
      k = 150, // scale
      x = 480, y = 250, // translate
      dx, dy, lambda = 0, phi = 0, // center
      deltaLambda = 0, deltaPhi = 0, deltaGamma = 0, rotate, projectRotate, // rotate
      theta = null, preclip = clip_antimeridian, // clip angle
      x0 = null, y0, x1, y1, postclip = identity, // clip extent
      delta2 = 0.5, projectResample = resample(projectTransform, delta2), // precision
      cache,
      cacheStream;

  function projection(point) {
    point = projectRotate(point[0] * radians, point[1] * radians);
    return [point[0] * k + dx, dy - point[1] * k];
  }

  function invert(point) {
    point = projectRotate.invert((point[0] - dx) / k, (dy - point[1]) / k);
    return point && [point[0] * degrees, point[1] * degrees];
  }

  function projectTransform(x, y) {
    return x = project(x, y), [x[0] * k + dx, dy - x[1] * k];
  }

  projection.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = transformRadians(preclip(rotate, projectResample(postclip(cacheStream = stream))));
  };

  projection.clipAngle = function(_) {
    return arguments.length ? (preclip = +_ ? clip_circle(theta = _ * radians, 6 * radians) : (theta = null, clip_antimeridian), reset()) : theta * degrees;
  };

  projection.clipExtent = function(_) {
    return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, identity) : extent_clipExtent(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  projection.scale = function(_) {
    return arguments.length ? (k = +_, recenter()) : k;
  };

  projection.translate = function(_) {
    return arguments.length ? (x = +_[0], y = +_[1], recenter()) : [x, y];
  };

  projection.center = function(_) {
    return arguments.length ? (lambda = _[0] % 360 * radians, phi = _[1] % 360 * radians, recenter()) : [lambda * degrees, phi * degrees];
  };

  projection.rotate = function(_) {
    return arguments.length ? (deltaLambda = _[0] % 360 * radians, deltaPhi = _[1] % 360 * radians, deltaGamma = _.length > 2 ? _[2] % 360 * radians : 0, recenter()) : [deltaLambda * degrees, deltaPhi * degrees, deltaGamma * degrees];
  };

  projection.precision = function(_) {
    return arguments.length ? (projectResample = resample(projectTransform, delta2 = _ * _), reset()) : sqrt(delta2);
  };

  projection.fitExtent = function(extent, object) {
    return fitExtent(projection, extent, object);
  };

  projection.fitSize = function(size, object) {
    return fitSize(projection, size, object);
  };

  function recenter() {
    projectRotate = compose(rotate = rotateRadians(deltaLambda, deltaPhi, deltaGamma), project);
    var center = project(lambda, phi);
    dx = x - center[0] * k;
    dy = y + center[1] * k;
    return reset();
  }

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return function() {
    project = projectAt.apply(this, arguments);
    projection.invert = project.invert && invert;
    return recenter();
  };
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/conic.js



function conicProjection(projectAt) {
  var phi0 = 0,
      phi1 = pi / 3,
      m = projectionMutator(projectAt),
      p = m(phi0, phi1);

  p.parallels = function(_) {
    return arguments.length ? m(phi0 = _[0] * radians, phi1 = _[1] * radians) : [phi0 * degrees, phi1 * degrees];
  };

  return p;
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/cylindricalEqualArea.js


function cylindricalEqualAreaRaw(phi0) {
  var cosPhi0 = cos(phi0);

  function forward(lambda, phi) {
    return [lambda * cosPhi0, sin(phi) / cosPhi0];
  }

  forward.invert = function(x, y) {
    return [x / cosPhi0, asin(y * cosPhi0)];
  };

  return forward;
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/conicEqualArea.js




function conicEqualAreaRaw(y0, y1) {
  var sy0 = sin(y0), n = (sy0 + sin(y1)) / 2;

  // Are the parallels symmetrical around the Equator?
  if (abs(n) < epsilon) return cylindricalEqualAreaRaw(y0);

  var c = 1 + sy0 * (2 * n - sy0), r0 = sqrt(c) / n;

  function project(x, y) {
    var r = sqrt(c - 2 * n * sin(y)) / n;
    return [r * sin(x *= n), r0 - r * cos(x)];
  }

  project.invert = function(x, y) {
    var r0y = r0 - y;
    return [atan2(x, abs(r0y)) / n * math_sign(r0y), asin((c - (x * x + r0y * r0y) * n * n) / (2 * n))];
  };

  return project;
}

/* harmony default export */ var conicEqualArea = (function() {
  return conicProjection(conicEqualAreaRaw)
      .scale(155.424)
      .center([0, 33.6442]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/albers.js


/* harmony default export */ var albers = (function() {
  return conicEqualArea()
      .parallels([29.5, 45.5])
      .scale(1070)
      .translate([480, 250])
      .rotate([96, 0])
      .center([-0.6, 38.7]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/albersUsa.js





// The projections must have mutually exclusive clip regions on the sphere,
// as this will avoid emitting interleaving lines and polygons.
function multiplex(streams) {
  var n = streams.length;
  return {
    point: function(x, y) { var i = -1; while (++i < n) streams[i].point(x, y); },
    sphere: function() { var i = -1; while (++i < n) streams[i].sphere(); },
    lineStart: function() { var i = -1; while (++i < n) streams[i].lineStart(); },
    lineEnd: function() { var i = -1; while (++i < n) streams[i].lineEnd(); },
    polygonStart: function() { var i = -1; while (++i < n) streams[i].polygonStart(); },
    polygonEnd: function() { var i = -1; while (++i < n) streams[i].polygonEnd(); }
  };
}

// A composite projection for the United States, configured by default for
// 960×500. The projection also works quite well at 960×600 if you change the
// scale to 1285 and adjust the translate accordingly. The set of standard
// parallels for each region comes from USGS, which is published here:
// http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html#albers
/* harmony default export */ var projection_albersUsa = (function() {
  var cache,
      cacheStream,
      lower48 = albers(), lower48Point,
      alaska = conicEqualArea().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]), alaskaPoint, // EPSG:3338
      hawaii = conicEqualArea().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]), hawaiiPoint, // ESRI:102007
      point, pointStream = {point: function(x, y) { point = [x, y]; }};

  function albersUsa(coordinates) {
    var x = coordinates[0], y = coordinates[1];
    return point = null,
        (lower48Point.point(x, y), point)
        || (alaskaPoint.point(x, y), point)
        || (hawaiiPoint.point(x, y), point);
  }

  albersUsa.invert = function(coordinates) {
    var k = lower48.scale(),
        t = lower48.translate(),
        x = (coordinates[0] - t[0]) / k,
        y = (coordinates[1] - t[1]) / k;
    return (y >= 0.120 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska
        : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii
        : lower48).invert(coordinates);
  };

  albersUsa.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex([lower48.stream(cacheStream = stream), alaska.stream(stream), hawaii.stream(stream)]);
  };

  albersUsa.precision = function(_) {
    if (!arguments.length) return lower48.precision();
    lower48.precision(_), alaska.precision(_), hawaii.precision(_);
    return reset();
  };

  albersUsa.scale = function(_) {
    if (!arguments.length) return lower48.scale();
    lower48.scale(_), alaska.scale(_ * 0.35), hawaii.scale(_);
    return albersUsa.translate(lower48.translate());
  };

  albersUsa.translate = function(_) {
    if (!arguments.length) return lower48.translate();
    var k = lower48.scale(), x = +_[0], y = +_[1];

    lower48Point = lower48
        .translate(_)
        .clipExtent([[x - 0.455 * k, y - 0.238 * k], [x + 0.455 * k, y + 0.238 * k]])
        .stream(pointStream);

    alaskaPoint = alaska
        .translate([x - 0.307 * k, y + 0.201 * k])
        .clipExtent([[x - 0.425 * k + epsilon, y + 0.120 * k + epsilon], [x - 0.214 * k - epsilon, y + 0.234 * k - epsilon]])
        .stream(pointStream);

    hawaiiPoint = hawaii
        .translate([x - 0.205 * k, y + 0.212 * k])
        .clipExtent([[x - 0.214 * k + epsilon, y + 0.166 * k + epsilon], [x - 0.115 * k - epsilon, y + 0.234 * k - epsilon]])
        .stream(pointStream);

    return reset();
  };

  albersUsa.fitExtent = function(extent, object) {
    return fitExtent(albersUsa, extent, object);
  };

  albersUsa.fitSize = function(size, object) {
    return fitSize(albersUsa, size, object);
  };

  function reset() {
    cache = cacheStream = null;
    return albersUsa;
  }

  return albersUsa.scale(1070);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/azimuthal.js


function azimuthalRaw(scale) {
  return function(x, y) {
    var cx = cos(x),
        cy = cos(y),
        k = scale(cx * cy);
    return [
      k * cy * sin(x),
      k * sin(y)
    ];
  }
}

function azimuthalInvert(angle) {
  return function(x, y) {
    var z = sqrt(x * x + y * y),
        c = angle(z),
        sc = sin(c),
        cc = cos(c);
    return [
      atan2(x * sc, z * cc),
      asin(z && y * sc / z)
    ];
  }
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/azimuthalEqualArea.js




var azimuthalEqualAreaRaw = azimuthalRaw(function(cxcy) {
  return sqrt(2 / (1 + cxcy));
});

azimuthalEqualAreaRaw.invert = azimuthalInvert(function(z) {
  return 2 * asin(z / 2);
});

/* harmony default export */ var azimuthalEqualArea = (function() {
  return projection_projection(azimuthalEqualAreaRaw)
      .scale(124.75)
      .clipAngle(180 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/azimuthalEquidistant.js




var azimuthalEquidistantRaw = azimuthalRaw(function(c) {
  return (c = acos(c)) && c / sin(c);
});

azimuthalEquidistantRaw.invert = azimuthalInvert(function(z) {
  return z;
});

/* harmony default export */ var azimuthalEquidistant = (function() {
  return projection_projection(azimuthalEquidistantRaw)
      .scale(79.4188)
      .clipAngle(180 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/mercator.js




function mercatorRaw(lambda, phi) {
  return [lambda, log(tan((halfPi + phi) / 2))];
}

mercatorRaw.invert = function(x, y) {
  return [x, 2 * atan(exp(y)) - halfPi];
};

/* harmony default export */ var mercator = (function() {
  return mercatorProjection(mercatorRaw)
      .scale(961 / tau);
});

function mercatorProjection(project) {
  var m = projection_projection(project),
      center = m.center,
      scale = m.scale,
      translate = m.translate,
      clipExtent = m.clipExtent,
      x0 = null, y0, x1, y1; // clip extent

  m.scale = function(_) {
    return arguments.length ? (scale(_), reclip()) : scale();
  };

  m.translate = function(_) {
    return arguments.length ? (translate(_), reclip()) : translate();
  };

  m.center = function(_) {
    return arguments.length ? (center(_), reclip()) : center();
  };

  m.clipExtent = function(_) {
    return arguments.length ? ((_ == null ? x0 = y0 = x1 = y1 = null : (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1])), reclip()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  function reclip() {
    var k = pi * scale(),
        t = m(src_rotation(m.rotate()).invert([0, 0]));
    return clipExtent(x0 == null
        ? [[t[0] - k, t[1] - k], [t[0] + k, t[1] + k]] : project === mercatorRaw
        ? [[Math.max(t[0] - k, x0), y0], [Math.min(t[0] + k, x1), y1]]
        : [[x0, Math.max(t[1] - k, y0)], [x1, Math.min(t[1] + k, y1)]]);
  }

  return reclip();
}

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/conicConformal.js




function tany(y) {
  return tan((halfPi + y) / 2);
}

function conicConformalRaw(y0, y1) {
  var cy0 = cos(y0),
      n = y0 === y1 ? sin(y0) : log(cy0 / cos(y1)) / log(tany(y1) / tany(y0)),
      f = cy0 * pow(tany(y0), n) / n;

  if (!n) return mercatorRaw;

  function project(x, y) {
    if (f > 0) { if (y < -halfPi + epsilon) y = -halfPi + epsilon; }
    else { if (y > halfPi - epsilon) y = halfPi - epsilon; }
    var r = f / pow(tany(y), n);
    return [r * sin(n * x), f - r * cos(n * x)];
  }

  project.invert = function(x, y) {
    var fy = f - y, r = math_sign(n) * sqrt(x * x + fy * fy);
    return [atan2(x, abs(fy)) / n * math_sign(fy), 2 * atan(pow(f / r, 1 / n)) - halfPi];
  };

  return project;
}

/* harmony default export */ var conicConformal = (function() {
  return conicProjection(conicConformalRaw)
      .scale(109.5)
      .parallels([30, 30]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/equirectangular.js


function equirectangularRaw(lambda, phi) {
  return [lambda, phi];
}

equirectangularRaw.invert = equirectangularRaw;

/* harmony default export */ var equirectangular = (function() {
  return projection_projection(equirectangularRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/conicEquidistant.js




function conicEquidistantRaw(y0, y1) {
  var cy0 = cos(y0),
      n = y0 === y1 ? sin(y0) : (cy0 - cos(y1)) / (y1 - y0),
      g = cy0 / n + y0;

  if (abs(n) < epsilon) return equirectangularRaw;

  function project(x, y) {
    var gy = g - y, nx = n * x;
    return [gy * sin(nx), g - gy * cos(nx)];
  }

  project.invert = function(x, y) {
    var gy = g - y;
    return [atan2(x, abs(gy)) / n * math_sign(gy), g - math_sign(n) * sqrt(x * x + gy * gy)];
  };

  return project;
}

/* harmony default export */ var conicEquidistant = (function() {
  return conicProjection(conicEquidistantRaw)
      .scale(131.154)
      .center([0, 13.9389]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/gnomonic.js




function gnomonicRaw(x, y) {
  var cy = cos(y), k = cos(x) * cy;
  return [cy * sin(x) / k, sin(y) / k];
}

gnomonicRaw.invert = azimuthalInvert(atan);

/* harmony default export */ var gnomonic = (function() {
  return projection_projection(gnomonicRaw)
      .scale(144.049)
      .clipAngle(60);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/identity.js





function scaleTranslate(kx, ky, tx, ty) {
  return kx === 1 && ky === 1 && tx === 0 && ty === 0 ? identity : transformer({
    point: function(x, y) {
      this.stream.point(x * kx + tx, y * ky + ty);
    }
  });
}

/* harmony default export */ var projection_identity = (function() {
  var k = 1, tx = 0, ty = 0, sx = 1, sy = 1, transform = identity, // scale, translate and reflect
      x0 = null, y0, x1, y1, clip = identity, // clip extent
      cache,
      cacheStream,
      projection;

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return projection = {
    stream: function(stream) {
      return cache && cacheStream === stream ? cache : cache = transform(clip(cacheStream = stream));
    },
    clipExtent: function(_) {
      return arguments.length ? (clip = _ == null ? (x0 = y0 = x1 = y1 = null, identity) : extent_clipExtent(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
    },
    scale: function(_) {
      return arguments.length ? (transform = scaleTranslate((k = +_) * sx, k * sy, tx, ty), reset()) : k;
    },
    translate: function(_) {
      return arguments.length ? (transform = scaleTranslate(k * sx, k * sy, tx = +_[0], ty = +_[1]), reset()) : [tx, ty];
    },
    reflectX: function(_) {
      return arguments.length ? (transform = scaleTranslate(k * (sx = _ ? -1 : 1), k * sy, tx, ty), reset()) : sx < 0;
    },
    reflectY: function(_) {
      return arguments.length ? (transform = scaleTranslate(k * sx, k * (sy = _ ? -1 : 1), tx, ty), reset()) : sy < 0;
    },
    fitExtent: function(extent, object) {
      return fitExtent(projection, extent, object);
    },
    fitSize: function(size, object) {
      return fitSize(projection, size, object);
    }
  };
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/orthographic.js




function orthographicRaw(x, y) {
  return [cos(y) * sin(x), sin(y)];
}

orthographicRaw.invert = azimuthalInvert(asin);

/* harmony default export */ var orthographic = (function() {
  return projection_projection(orthographicRaw)
      .scale(249.5)
      .clipAngle(90 + epsilon);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/stereographic.js




function stereographicRaw(x, y) {
  var cy = cos(y), k = 1 + cos(x) * cy;
  return [cy * sin(x) / k, sin(y) / k];
}

stereographicRaw.invert = azimuthalInvert(function(z) {
  return 2 * atan(z);
});

/* harmony default export */ var stereographic = (function() {
  return projection_projection(stereographicRaw)
      .scale(250)
      .clipAngle(142);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/src/projection/transverseMercator.js



function transverseMercatorRaw(lambda, phi) {
  return [log(tan((halfPi + phi) / 2)), -lambda];
}

transverseMercatorRaw.invert = function(x, y) {
  return [-y, 2 * atan(exp(x)) - halfPi];
};

/* harmony default export */ var transverseMercator = (function() {
  var m = mercatorProjection(transverseMercatorRaw),
      center = m.center,
      rotate = m.rotate;

  m.center = function(_) {
    return arguments.length ? center([-_[1], _[0]]) : (_ = center(), [_[1], -_[0]]);
  };

  m.rotate = function(_) {
    return arguments.length ? rotate([_[0], _[1], _.length > 2 ? _[2] + 90 : 90]) : (_ = rotate(), [_[0], _[1], _[2] - 90]);
  };

  return rotate([0, 0, 90])
      .scale(159.155);
});

// CONCATENATED MODULE: ./node_modules/d3-geo/index.js
/* concated harmony reexport geoArea */__webpack_require__.d(__webpack_exports__, "geoArea", function() { return src_area; });
/* concated harmony reexport geoBounds */__webpack_require__.d(__webpack_exports__, "geoBounds", function() { return bounds; });
/* concated harmony reexport geoCentroid */__webpack_require__.d(__webpack_exports__, "geoCentroid", function() { return centroid; });
/* concated harmony reexport geoCircle */__webpack_require__.d(__webpack_exports__, "geoCircle", function() { return src_circle; });
/* concated harmony reexport geoClipExtent */__webpack_require__.d(__webpack_exports__, "geoClipExtent", function() { return clip_extent; });
/* concated harmony reexport geoContains */__webpack_require__.d(__webpack_exports__, "geoContains", function() { return contains; });
/* concated harmony reexport geoDistance */__webpack_require__.d(__webpack_exports__, "geoDistance", function() { return distance; });
/* concated harmony reexport geoGraticule */__webpack_require__.d(__webpack_exports__, "geoGraticule", function() { return graticule_graticule; });
/* concated harmony reexport geoGraticule10 */__webpack_require__.d(__webpack_exports__, "geoGraticule10", function() { return graticule10; });
/* concated harmony reexport geoInterpolate */__webpack_require__.d(__webpack_exports__, "geoInterpolate", function() { return src_interpolate; });
/* concated harmony reexport geoLength */__webpack_require__.d(__webpack_exports__, "geoLength", function() { return src_length; });
/* concated harmony reexport geoPath */__webpack_require__.d(__webpack_exports__, "geoPath", function() { return src_path; });
/* concated harmony reexport geoAlbers */__webpack_require__.d(__webpack_exports__, "geoAlbers", function() { return albers; });
/* concated harmony reexport geoAlbersUsa */__webpack_require__.d(__webpack_exports__, "geoAlbersUsa", function() { return projection_albersUsa; });
/* concated harmony reexport geoAzimuthalEqualArea */__webpack_require__.d(__webpack_exports__, "geoAzimuthalEqualArea", function() { return azimuthalEqualArea; });
/* concated harmony reexport geoAzimuthalEqualAreaRaw */__webpack_require__.d(__webpack_exports__, "geoAzimuthalEqualAreaRaw", function() { return azimuthalEqualAreaRaw; });
/* concated harmony reexport geoAzimuthalEquidistant */__webpack_require__.d(__webpack_exports__, "geoAzimuthalEquidistant", function() { return azimuthalEquidistant; });
/* concated harmony reexport geoAzimuthalEquidistantRaw */__webpack_require__.d(__webpack_exports__, "geoAzimuthalEquidistantRaw", function() { return azimuthalEquidistantRaw; });
/* concated harmony reexport geoConicConformal */__webpack_require__.d(__webpack_exports__, "geoConicConformal", function() { return conicConformal; });
/* concated harmony reexport geoConicConformalRaw */__webpack_require__.d(__webpack_exports__, "geoConicConformalRaw", function() { return conicConformalRaw; });
/* concated harmony reexport geoConicEqualArea */__webpack_require__.d(__webpack_exports__, "geoConicEqualArea", function() { return conicEqualArea; });
/* concated harmony reexport geoConicEqualAreaRaw */__webpack_require__.d(__webpack_exports__, "geoConicEqualAreaRaw", function() { return conicEqualAreaRaw; });
/* concated harmony reexport geoConicEquidistant */__webpack_require__.d(__webpack_exports__, "geoConicEquidistant", function() { return conicEquidistant; });
/* concated harmony reexport geoConicEquidistantRaw */__webpack_require__.d(__webpack_exports__, "geoConicEquidistantRaw", function() { return conicEquidistantRaw; });
/* concated harmony reexport geoEquirectangular */__webpack_require__.d(__webpack_exports__, "geoEquirectangular", function() { return equirectangular; });
/* concated harmony reexport geoEquirectangularRaw */__webpack_require__.d(__webpack_exports__, "geoEquirectangularRaw", function() { return equirectangularRaw; });
/* concated harmony reexport geoGnomonic */__webpack_require__.d(__webpack_exports__, "geoGnomonic", function() { return gnomonic; });
/* concated harmony reexport geoGnomonicRaw */__webpack_require__.d(__webpack_exports__, "geoGnomonicRaw", function() { return gnomonicRaw; });
/* concated harmony reexport geoIdentity */__webpack_require__.d(__webpack_exports__, "geoIdentity", function() { return projection_identity; });
/* concated harmony reexport geoProjection */__webpack_require__.d(__webpack_exports__, "geoProjection", function() { return projection_projection; });
/* concated harmony reexport geoProjectionMutator */__webpack_require__.d(__webpack_exports__, "geoProjectionMutator", function() { return projectionMutator; });
/* concated harmony reexport geoMercator */__webpack_require__.d(__webpack_exports__, "geoMercator", function() { return mercator; });
/* concated harmony reexport geoMercatorRaw */__webpack_require__.d(__webpack_exports__, "geoMercatorRaw", function() { return mercatorRaw; });
/* concated harmony reexport geoOrthographic */__webpack_require__.d(__webpack_exports__, "geoOrthographic", function() { return orthographic; });
/* concated harmony reexport geoOrthographicRaw */__webpack_require__.d(__webpack_exports__, "geoOrthographicRaw", function() { return orthographicRaw; });
/* concated harmony reexport geoStereographic */__webpack_require__.d(__webpack_exports__, "geoStereographic", function() { return stereographic; });
/* concated harmony reexport geoStereographicRaw */__webpack_require__.d(__webpack_exports__, "geoStereographicRaw", function() { return stereographicRaw; });
/* concated harmony reexport geoTransverseMercator */__webpack_require__.d(__webpack_exports__, "geoTransverseMercator", function() { return transverseMercator; });
/* concated harmony reexport geoTransverseMercatorRaw */__webpack_require__.d(__webpack_exports__, "geoTransverseMercatorRaw", function() { return transverseMercatorRaw; });
/* concated harmony reexport geoRotation */__webpack_require__.d(__webpack_exports__, "geoRotation", function() { return src_rotation; });
/* concated harmony reexport geoStream */__webpack_require__.d(__webpack_exports__, "geoStream", function() { return src_stream; });
/* concated harmony reexport geoTransform */__webpack_require__.d(__webpack_exports__, "geoTransform", function() { return src_transform; });




 // DEPRECATED! Use d3.geoIdentity().clipExtent(…).


























/***/ }),
/* 1 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__1__;

/***/ }),
/* 2 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";

// CONCATENATED MODULE: ./node_modules/d3-array/src/ascending.js
/* harmony default export */ var ascending = (function(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/bisector.js


/* harmony default export */ var bisector = (function(compare) {
  if (compare.length === 1) compare = ascendingComparator(compare);
  return {
    left: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) < 0) lo = mid + 1;
        else hi = mid;
      }
      return lo;
    },
    right: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) > 0) hi = mid;
        else lo = mid + 1;
      }
      return lo;
    }
  };
});

function ascendingComparator(f) {
  return function(d, x) {
    return ascending(f(d), x);
  };
}

// CONCATENATED MODULE: ./node_modules/d3-array/src/bisect.js



var ascendingBisect = bisector(ascending);
var bisectRight = ascendingBisect.right;
var bisectLeft = ascendingBisect.left;
/* harmony default export */ var bisect = (bisectRight);

// CONCATENATED MODULE: ./node_modules/d3-array/src/pairs.js
/* harmony default export */ var pairs = (function(array, f) {
  if (f == null) f = pair;
  var i = 0, n = array.length - 1, p = array[0], pairs = new Array(n < 0 ? 0 : n);
  while (i < n) pairs[i] = f(p, p = array[++i]);
  return pairs;
});

function pair(a, b) {
  return [a, b];
}

// CONCATENATED MODULE: ./node_modules/d3-array/src/cross.js


/* harmony default export */ var cross = (function(values0, values1, reduce) {
  var n0 = values0.length,
      n1 = values1.length,
      values = new Array(n0 * n1),
      i0,
      i1,
      i,
      value0;

  if (reduce == null) reduce = pair;

  for (i0 = i = 0; i0 < n0; ++i0) {
    for (value0 = values0[i0], i1 = 0; i1 < n1; ++i1, ++i) {
      values[i] = reduce(value0, values1[i1]);
    }
  }

  return values;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/descending.js
/* harmony default export */ var descending = (function(a, b) {
  return b < a ? -1 : b > a ? 1 : b >= a ? 0 : NaN;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/number.js
/* harmony default export */ var number = (function(x) {
  return x === null ? NaN : +x;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/variance.js


/* harmony default export */ var variance = (function(values, valueof) {
  var n = values.length,
      m = 0,
      i = -1,
      mean = 0,
      value,
      delta,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  if (m > 1) return sum / (m - 1);
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/deviation.js


/* harmony default export */ var deviation = (function(array, f) {
  var v = variance(array, f);
  return v ? Math.sqrt(v) : v;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/extent.js
/* harmony default export */ var extent = (function(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  return [min, max];
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/array.js
var array_array = Array.prototype;

var slice = array_array.slice;
var map = array_array.map;

// CONCATENATED MODULE: ./node_modules/d3-array/src/constant.js
/* harmony default export */ var constant = (function(x) {
  return function() {
    return x;
  };
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/identity.js
/* harmony default export */ var identity = (function(x) {
  return x;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/range.js
/* harmony default export */ var range = (function(start, stop, step) {
  start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;

  var i = -1,
      n = Math.max(0, Math.ceil((stop - start) / step)) | 0,
      range = new Array(n);

  while (++i < n) {
    range[i] = start + i * step;
  }

  return range;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/ticks.js
var e10 = Math.sqrt(50),
    e5 = Math.sqrt(10),
    e2 = Math.sqrt(2);

/* harmony default export */ var ticks = (function(start, stop, count) {
  var reverse,
      i = -1,
      n,
      ticks,
      step;

  stop = +stop, start = +start, count = +count;
  if (start === stop && count > 0) return [start];
  if (reverse = stop < start) n = start, start = stop, stop = n;
  if ((step = tickIncrement(start, stop, count)) === 0 || !isFinite(step)) return [];

  if (step > 0) {
    start = Math.ceil(start / step);
    stop = Math.floor(stop / step);
    ticks = new Array(n = Math.ceil(stop - start + 1));
    while (++i < n) ticks[i] = (start + i) * step;
  } else {
    start = Math.floor(start * step);
    stop = Math.ceil(stop * step);
    ticks = new Array(n = Math.ceil(start - stop + 1));
    while (++i < n) ticks[i] = (start - i) / step;
  }

  if (reverse) ticks.reverse();

  return ticks;
});

function tickIncrement(start, stop, count) {
  var step = (stop - start) / Math.max(0, count),
      power = Math.floor(Math.log(step) / Math.LN10),
      error = step / Math.pow(10, power);
  return power >= 0
      ? (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1) * Math.pow(10, power)
      : -Math.pow(10, -power) / (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1);
}

function tickStep(start, stop, count) {
  var step0 = Math.abs(stop - start) / Math.max(0, count),
      step1 = Math.pow(10, Math.floor(Math.log(step0) / Math.LN10)),
      error = step0 / step1;
  if (error >= e10) step1 *= 10;
  else if (error >= e5) step1 *= 5;
  else if (error >= e2) step1 *= 2;
  return stop < start ? -step1 : step1;
}

// CONCATENATED MODULE: ./node_modules/d3-array/src/threshold/sturges.js
/* harmony default export */ var sturges = (function(values) {
  return Math.ceil(Math.log(values.length) / Math.LN2) + 1;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/histogram.js









/* harmony default export */ var src_histogram = (function() {
  var value = identity,
      domain = extent,
      threshold = sturges;

  function histogram(data) {
    var i,
        n = data.length,
        x,
        values = new Array(n);

    for (i = 0; i < n; ++i) {
      values[i] = value(data[i], i, data);
    }

    var xz = domain(values),
        x0 = xz[0],
        x1 = xz[1],
        tz = threshold(values, x0, x1);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      tz = tickStep(x0, x1, tz);
      tz = range(Math.ceil(x0 / tz) * tz, x1, tz); // exclusive
    }

    // Remove any thresholds outside the domain.
    var m = tz.length;
    while (tz[0] <= x0) tz.shift(), --m;
    while (tz[m - 1] > x1) tz.pop(), --m;

    var bins = new Array(m + 1),
        bin;

    // Initialize bins.
    for (i = 0; i <= m; ++i) {
      bin = bins[i] = [];
      bin.x0 = i > 0 ? tz[i - 1] : x0;
      bin.x1 = i < m ? tz[i] : x1;
    }

    // Assign data to bins by value, ignoring any outside the domain.
    for (i = 0; i < n; ++i) {
      x = values[i];
      if (x0 <= x && x <= x1) {
        bins[bisect(tz, x, 0, m)].push(data[i]);
      }
    }

    return bins;
  }

  histogram.value = function(_) {
    return arguments.length ? (value = typeof _ === "function" ? _ : constant(_), histogram) : value;
  };

  histogram.domain = function(_) {
    return arguments.length ? (domain = typeof _ === "function" ? _ : constant([_[0], _[1]]), histogram) : domain;
  };

  histogram.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), histogram) : threshold;
  };

  return histogram;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/quantile.js


/* harmony default export */ var quantile = (function(values, p, valueof) {
  if (valueof == null) valueof = number;
  if (!(n = values.length)) return;
  if ((p = +p) <= 0 || n < 2) return +valueof(values[0], 0, values);
  if (p >= 1) return +valueof(values[n - 1], n - 1, values);
  var n,
      i = (n - 1) * p,
      i0 = Math.floor(i),
      value0 = +valueof(values[i0], i0, values),
      value1 = +valueof(values[i0 + 1], i0 + 1, values);
  return value0 + (value1 - value0) * (i - i0);
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/threshold/freedmanDiaconis.js





/* harmony default export */ var freedmanDiaconis = (function(values, min, max) {
  values = map.call(values, number).sort(ascending);
  return Math.ceil((max - min) / (2 * (quantile(values, 0.75) - quantile(values, 0.25)) * Math.pow(values.length, -1 / 3)));
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/threshold/scott.js


/* harmony default export */ var scott = (function(values, min, max) {
  return Math.ceil((max - min) / (3.5 * deviation(values) * Math.pow(values.length, -1 / 3)));
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/max.js
/* harmony default export */ var src_max = (function(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  return max;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/mean.js


/* harmony default export */ var src_mean = (function(values, valueof) {
  var n = values.length,
      m = n,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) sum += value;
      else --m;
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) sum += value;
      else --m;
    }
  }

  if (m) return sum / m;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/median.js




/* harmony default export */ var median = (function(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      numbers = [];

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) {
        numbers.push(value);
      }
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) {
        numbers.push(value);
      }
    }
  }

  return quantile(numbers.sort(ascending), 0.5);
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/merge.js
/* harmony default export */ var merge = (function(arrays) {
  var n = arrays.length,
      m,
      i = -1,
      j = 0,
      merged,
      array;

  while (++i < n) j += arrays[i].length;
  merged = new Array(j);

  while (--n >= 0) {
    array = arrays[n];
    m = array.length;
    while (--m >= 0) {
      merged[--j] = array[m];
    }
  }

  return merged;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/min.js
/* harmony default export */ var src_min = (function(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  return min;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/permute.js
/* harmony default export */ var permute = (function(array, indexes) {
  var i = indexes.length, permutes = new Array(i);
  while (i--) permutes[i] = array[indexes[i]];
  return permutes;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/scan.js


/* harmony default export */ var scan = (function(values, compare) {
  if (!(n = values.length)) return;
  var n,
      i = 0,
      j = 0,
      xi,
      xj = values[j];

  if (compare == null) compare = ascending;

  while (++i < n) {
    if (compare(xi = values[i], xj) < 0 || compare(xj, xj) !== 0) {
      xj = xi, j = i;
    }
  }

  if (compare(xj, xj) === 0) return j;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/shuffle.js
/* harmony default export */ var shuffle = (function(array, i0, i1) {
  var m = (i1 == null ? array.length : i1) - (i0 = i0 == null ? 0 : +i0),
      t,
      i;

  while (m) {
    i = Math.random() * m-- | 0;
    t = array[m + i0];
    array[m + i0] = array[i + i0];
    array[i + i0] = t;
  }

  return array;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/sum.js
/* harmony default export */ var src_sum = (function(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (value = +values[i]) sum += value; // Note: zero and null are equivalent.
    }
  }

  else {
    while (++i < n) {
      if (value = +valueof(values[i], i, values)) sum += value;
    }
  }

  return sum;
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/transpose.js


/* harmony default export */ var src_transpose = (function(matrix) {
  if (!(n = matrix.length)) return [];
  for (var i = -1, m = src_min(matrix, transpose_length), transpose = new Array(m); ++i < m;) {
    for (var j = -1, n, row = transpose[i] = new Array(n); ++j < n;) {
      row[j] = matrix[j][i];
    }
  }
  return transpose;
});

function transpose_length(d) {
  return d.length;
}

// CONCATENATED MODULE: ./node_modules/d3-array/src/zip.js


/* harmony default export */ var zip = (function() {
  return src_transpose(arguments);
});

// CONCATENATED MODULE: ./node_modules/d3-array/src/index.js
/* unused concated harmony import bisect */
/* unused concated harmony import bisectRight */
/* unused concated harmony import bisectLeft */
/* unused concated harmony import ascending */
/* unused concated harmony import bisector */
/* unused concated harmony import cross */
/* unused concated harmony import descending */
/* unused concated harmony import deviation */
/* unused concated harmony import extent */
/* unused concated harmony import histogram */
/* unused concated harmony import thresholdFreedmanDiaconis */
/* unused concated harmony import thresholdScott */
/* unused concated harmony import thresholdSturges */
/* unused concated harmony import max */
/* unused concated harmony import mean */
/* unused concated harmony import median */
/* concated harmony reexport merge */__webpack_require__.d(__webpack_exports__, "a", function() { return merge; });
/* unused concated harmony import min */
/* unused concated harmony import pairs */
/* unused concated harmony import permute */
/* unused concated harmony import quantile */
/* concated harmony reexport range */__webpack_require__.d(__webpack_exports__, "b", function() { return range; });
/* unused concated harmony import scan */
/* unused concated harmony import shuffle */
/* unused concated harmony import sum */
/* unused concated harmony import ticks */
/* unused concated harmony import tickIncrement */
/* unused concated harmony import tickStep */
/* unused concated harmony import transpose */
/* unused concated harmony import variance */
/* unused concated harmony import zip */





























/***/ }),
/* 3 */
/***/ (function(module, exports) {

module.exports = __WEBPACK_EXTERNAL_MODULE__3__;

/***/ }),
/* 4 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _ComposableMap = __webpack_require__(8);

Object.defineProperty(exports, "ComposableMap", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_ComposableMap).default;
  }
});

var _ZoomableGroup = __webpack_require__(10);

Object.defineProperty(exports, "ZoomableGroup", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_ZoomableGroup).default;
  }
});

var _ZoomableGlobe = __webpack_require__(11);

Object.defineProperty(exports, "ZoomableGlobe", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_ZoomableGlobe).default;
  }
});

var _Geographies = __webpack_require__(12);

Object.defineProperty(exports, "Geographies", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Geographies).default;
  }
});

var _Geography = __webpack_require__(13);

Object.defineProperty(exports, "Geography", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Geography).default;
  }
});

var _Marker = __webpack_require__(14);

Object.defineProperty(exports, "Marker", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Marker).default;
  }
});

var _Markers = __webpack_require__(15);

Object.defineProperty(exports, "Markers", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Markers).default;
  }
});

var _Line = __webpack_require__(16);

Object.defineProperty(exports, "Line", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Line).default;
  }
});

var _Lines = __webpack_require__(17);

Object.defineProperty(exports, "Lines", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Lines).default;
  }
});

var _MapGroup = __webpack_require__(6);

Object.defineProperty(exports, "MapGroup", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_MapGroup).default;
  }
});

var _Annotation = __webpack_require__(18);

Object.defineProperty(exports, "Annotation", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Annotation).default;
  }
});

var _Annotations = __webpack_require__(19);

Object.defineProperty(exports, "Annotations", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Annotations).default;
  }
});

var _Graticule = __webpack_require__(20);

Object.defineProperty(exports, "Graticule", {
  enumerable: true,
  get: function get() {
    return _interopRequireDefault(_Graticule).default;
  }
});

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

/***/ }),
/* 5 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.calculateResizeFactor = calculateResizeFactor;
exports.calculateMousePosition = calculateMousePosition;
exports.isChildOfType = isChildOfType;
exports.createNewChildren = createNewChildren;
exports.roundPath = roundPath;
exports.createConnectorPath = createConnectorPath;
exports.createTextAnchor = createTextAnchor;
exports.computeBackdrop = computeBackdrop;

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function calculateResizeFactor(actualDimension, baseDimension) {
  if (actualDimension === 0) return 1;
  return 1 / 100 * (100 / actualDimension * baseDimension);
}

function calculateMousePosition(direction, projection, props, zoom, resizeFactor) {
  var center = arguments.length > 5 && arguments[5] !== undefined ? arguments[5] : props.center;
  var width = arguments.length > 6 && arguments[6] !== undefined ? arguments[6] : props.width;
  var height = arguments.length > 7 && arguments[7] !== undefined ? arguments[7] : props.height;

  var reference = { x: 0, y: 1 };
  var canRotate = !!projection.rotate;
  var reverseRotation = !!canRotate ? projection.rotate().map(function (item) {
    return -item;
  }) : false;
  var point = !!reverseRotation ? projection.rotate(reverseRotation)([-center[0], -center[1]]) : projection([center[0], center[1]]);
  var returner = point ? (point[reference[direction]] - (reference[direction] === 0 ? width : height) / 2) * zoom * (1 / resizeFactor) : 0;
  if (canRotate) projection.rotate([-reverseRotation[0], -reverseRotation[1], -reverseRotation[2]]);
  return !!reverseRotation ? returner : -returner;
}

function isChildOfType(child, expectedTypes) {
  return expectedTypes.indexOf(child.props.componentIdentifier) !== -1;
}

function createNewChildren(children, props) {
  if (!children) return;
  if (!children.length) {
    return isChildOfType(children, ["Geographies"]) ? _react2.default.cloneElement(children, {
      projection: props.projection
    }) : isChildOfType(children, ["Group", "Markers", "Lines", "Annotations", "Annotation", "Graticule"]) ? _react2.default.cloneElement(children, {
      projection: props.projection,
      zoom: props.zoom,
      width: props.width,
      height: props.height,
      groupName: props.groupName,
      itemName: props.itemName
    }) : children;
  } else {
    return children.map(function (child, i) {
      if (!child) return;
      return isChildOfType(child, ["Geographies"]) ? _react2.default.cloneElement(child, {
        key: "zoomable-child-" + i,
        projection: props.projection
      }) : isChildOfType(child, ["Group", "Markers", "Lines", "Annotations", "Annotation", "Graticule"]) ? _react2.default.cloneElement(child, {
        key: "zoomable-child-" + i,
        projection: props.projection,
        zoom: props.zoom,
        width: props.width,
        height: props.height,
        groupName: props.groupName,
        itemName: props.itemName
      }) : child;
    });
  }
}

function roundPath(path, precision) {
  if (!path) return;
  var query = /[\d\.-][\d\.e-]*/g;
  return path.replace(query, function (n) {
    return Math.round(n * (1 / precision)) / (1 / precision);
  });
}

function createConnectorPath(connectorType, endPoint, curve) {
  var e0 = endPoint[0];
  var e1 = endPoint[1];
  return "M0,0 Q " + (curve + 1) / 2 * e0 + "," + (e1 - (curve + 1) / 2 * e1) + " " + e0 + "," + e1;
}

function createTextAnchor(dx) {
  if (dx > 0) return "start";else if (dx < 0) return "end";else return "middle";
}

function computeBackdrop(projection, backdrop) {
  var canRotate = projection.rotate;
  var originalRotation = canRotate ? projection.rotate() : null;

  var tl = canRotate ? projection.rotate([0, 0, 0])([backdrop.x[0], backdrop.y[0]]) : projection([backdrop.x[0], backdrop.y[0]]);

  var br = canRotate ? projection.rotate([0, 0, 0])([backdrop.x[1], backdrop.y[1]]) : projection([backdrop.x[1], backdrop.y[1]]);

  var x = tl ? tl[0] : 0;
  var x0 = br ? br[0] : 0;

  var y = tl ? tl[1] : 0;
  var y0 = br ? br[1] : 0;

  var width = x0 - x;
  var height = y0 - y;

  if (originalRotation) projection.rotate(originalRotation);

  return { x: x, y: y, width: width, height: height };
}

/***/ }),
/* 6 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var MapGroup = function (_Component) {
  _inherits(MapGroup, _Component);

  function MapGroup() {
    _classCallCheck(this, MapGroup);

    return _possibleConstructorReturn(this, (MapGroup.__proto__ || Object.getPrototypeOf(MapGroup)).apply(this, arguments));
  }

  _createClass(MapGroup, [{
    key: "render",
    value: function render() {
      var _props = this.props,
          children = _props.children,
          projection = _props.projection,
          style = _props.style,
          zoom = _props.zoom,
          width = _props.width,
          height = _props.height,
          groupName = _props.groupName,
          itemName = _props.itemName;

      return _react2.default.createElement(
        "g",
        { className: "rsm-" + groupName, style: style },
        !children ? null : children.length === undefined ? _react2.default.cloneElement(children, {
          projection: projection,
          zoom: zoom,
          width: width,
          height: height
        }) : children.map(function (child, i) {
          return !child ? null : _react2.default.cloneElement(child, {
            key: child.key || itemName + "-" + i,
            projection: projection,
            zoom: zoom,
            width: width,
            height: height
          });
        })
      );
    }
  }]);

  return MapGroup;
}(_react.Component);

MapGroup.defaultProps = {
  componentIdentifier: "Group",
  groupName: "group",
  itemName: "group-item"
};

exports.default = MapGroup;

/***/ }),
/* 7 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.default = {
  scale: 160,
  xOffset: 0,
  yOffset: 0,
  rotation: [0, 0, 0],
  precision: 0.1
};

/***/ }),
/* 8 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _projections = __webpack_require__(9);

var _projections2 = _interopRequireDefault(_projections);

var _projectionConfig = __webpack_require__(7);

var _projectionConfig2 = _interopRequireDefault(_projectionConfig);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var ComposableMap = function (_Component) {
  _inherits(ComposableMap, _Component);

  function ComposableMap() {
    _classCallCheck(this, ComposableMap);

    var _this = _possibleConstructorReturn(this, (ComposableMap.__proto__ || Object.getPrototypeOf(ComposableMap)).call(this));

    _this.projection = _this.projection.bind(_this);
    return _this;
  }

  _createClass(ComposableMap, [{
    key: "projection",
    value: function projection() {
      var _props = this.props,
          projection = _props.projection,
          projectionConfig = _props.projectionConfig,
          width = _props.width,
          height = _props.height;


      return typeof projection !== "function" ? (0, _projections2.default)(width, height, projectionConfig, projection) : projection(width, height, projectionConfig);
    }
  }, {
    key: "render",
    value: function render() {
      var _props2 = this.props,
          width = _props2.width,
          height = _props2.height,
          style = _props2.style,
          className = _props2.className,
          showCenter = _props2.showCenter,
          children = _props2.children,
          aspectRatio = _props2.aspectRatio,
          viewBox = _props2.viewBox,
          defs = _props2.defs;


      return _react2.default.createElement(
        "svg",
        { width: width,
          height: height,
          viewBox: viewBox ? viewBox : "0 0 " + width + " " + height,
          className: "rsm-svg " + (className || ''),
          style: style,
          preserveAspectRatio: aspectRatio },
        defs && _react2.default.createElement(
          "defs",
          null,
          defs
        ),
        _react2.default.cloneElement(this.props.children, {
          projection: this.projection(),
          width: width,
          height: height
        }),
        showCenter && _react2.default.createElement(
          "g",
          null,
          _react2.default.createElement("rect", { x: width / 2 - 0.5, y: 0, width: 1, height: height, fill: "#e91e63" }),
          _react2.default.createElement("rect", { x: 0, y: height / 2 - 0.5, width: width, height: 1, fill: "#e91e63" })
        )
      );
    }
  }]);

  return ComposableMap;
}(_react.Component);

ComposableMap.defaultProps = {
  width: 800,
  height: 450,
  projection: "times",
  projectionConfig: _projectionConfig2.default,
  aspectRatio: "xMidYMid",
  viewBox: null
};

exports.default = ComposableMap;

/***/ }),
/* 9 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

exports.default = function (width, height, config, projectionName) {

  var scale = config.scale || _projectionConfig2.default.scale;
  var xOffset = config.xOffset || _projectionConfig2.default.xOffset;
  var yOffset = config.yOffset || _projectionConfig2.default.yOffset;
  var rotation = config.rotation || _projectionConfig2.default.rotation;
  var precision = config.precision || _projectionConfig2.default.precision;

  var baseProjection = projectionReference[projectionName]().scale(scale).translate([xOffset + width / 2, yOffset + height / 2]).precision(precision);

  return projectionName === "albersUsa" ? baseProjection : projectionName === "orthographic" ? baseProjection.rotate(rotation).clipAngle(90) : baseProjection.rotate(rotation);
};

var _d3GeoProjection = __webpack_require__(21);

var _d3Geo = __webpack_require__(0);

var _projectionConfig = __webpack_require__(7);

var _projectionConfig2 = _interopRequireDefault(_projectionConfig);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

var projectionReference = {
  mercator: _d3Geo.geoMercator,
  miller: _d3GeoProjection.geoMiller,
  times: _d3GeoProjection.geoTimes,
  robinson: _d3GeoProjection.geoRobinson,
  winkel3: _d3GeoProjection.geoWinkel3,
  eckert4: _d3GeoProjection.geoEckert4,
  albersUsa: _d3Geo.geoAlbersUsa,
  orthographic: _d3Geo.geoOrthographic
};

/***/ }),
/* 10 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

var _utils = __webpack_require__(5);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var ZoomableGroup = function (_Component) {
  _inherits(ZoomableGroup, _Component);

  function ZoomableGroup(props) {
    _classCallCheck(this, ZoomableGroup);

    var _this = _possibleConstructorReturn(this, (ZoomableGroup.__proto__ || Object.getPrototypeOf(ZoomableGroup)).call(this, props));

    var backdrop = (0, _utils.computeBackdrop)(props.projection, props.backdrop);

    _this.state = {
      mouseX: (0, _utils.calculateMousePosition)("x", props.projection, props, props.zoom, 1),
      mouseY: (0, _utils.calculateMousePosition)("y", props.projection, props, props.zoom, 1),
      mouseXStart: 0,
      mouseYStart: 0,
      isPressed: false,
      resizeFactorX: 1,
      resizeFactorY: 1,
      backdrop: {
        width: Math.round(backdrop.width),
        height: Math.round(backdrop.height),
        x: Math.round(backdrop.x),
        y: Math.round(backdrop.y)
      }
    };

    _this.handleMouseMove = _this.handleMouseMove.bind(_this);
    _this.handleMouseUp = _this.handleMouseUp.bind(_this);
    _this.handleMouseDown = _this.handleMouseDown.bind(_this);
    _this.handleTouchStart = _this.handleTouchStart.bind(_this);
    _this.handleTouchMove = _this.handleTouchMove.bind(_this);
    _this.handleResize = _this.handleResize.bind(_this);

    return _this;
  }

  _createClass(ZoomableGroup, [{
    key: "handleMouseMove",
    value: function handleMouseMove(_ref) {
      var pageX = _ref.pageX,
          pageY = _ref.pageY;

      if (this.props.disablePanning) return;
      if (!this.state.isPressed) return;
      this.setState({
        mouseX: pageX - this.state.mouseXStart,
        mouseY: pageY - this.state.mouseYStart
      });
    }
  }, {
    key: "handleTouchMove",
    value: function handleTouchMove(_ref2) {
      var touches = _ref2.touches;

      this.handleMouseMove(touches[0]);
    }
  }, {
    key: "handleMouseUp",
    value: function handleMouseUp() {
      if (this.props.disablePanning) return;
      if (!this.state.isPressed) return;
      this.setState({
        isPressed: false
      });
      if (!this.props.onMoveEnd) return;
      var _state = this.state,
          mouseX = _state.mouseX,
          mouseY = _state.mouseY,
          resizeFactorX = _state.resizeFactorX,
          resizeFactorY = _state.resizeFactorY;
      var _props = this.props,
          zoom = _props.zoom,
          width = _props.width,
          height = _props.height,
          projection = _props.projection,
          onMoveEnd = _props.onMoveEnd;

      var x = width / 2 - mouseX * resizeFactorX / zoom;
      var y = height / 2 - mouseY * resizeFactorY / zoom;
      var newCenter = projection.invert([x, y]);
      onMoveEnd(newCenter);
    }
  }, {
    key: "handleMouseDown",
    value: function handleMouseDown(_ref3) {
      var pageX = _ref3.pageX,
          pageY = _ref3.pageY;

      if (this.props.disablePanning) return;
      var _state2 = this.state,
          mouseX = _state2.mouseX,
          mouseY = _state2.mouseY,
          resizeFactorX = _state2.resizeFactorX,
          resizeFactorY = _state2.resizeFactorY;
      var _props2 = this.props,
          zoom = _props2.zoom,
          width = _props2.width,
          height = _props2.height,
          projection = _props2.projection,
          onMoveStart = _props2.onMoveStart;

      this.setState({
        isPressed: true,
        mouseXStart: pageX - mouseX,
        mouseYStart: pageY - mouseY
      });
      if (!onMoveStart) return;
      var x = width / 2 - mouseX * resizeFactorX / zoom;
      var y = height / 2 - mouseY * resizeFactorY / zoom;
      var currentCenter = projection.invert([x, y]);
      onMoveStart(currentCenter);
    }
  }, {
    key: "handleTouchStart",
    value: function handleTouchStart(_ref4) {
      var touches = _ref4.touches;

      if (touches.length > 1) {
        this.handleMouseDown(touches[0]);
      } else {
        this.handleMouseUp();
      }
    }
  }, {
    key: "preventTouchScroll",
    value: function preventTouchScroll(evt) {
      if (evt.touches.length > 1) {
        evt.preventDefault();
      }
    }
  }, {
    key: "componentWillReceiveProps",
    value: function componentWillReceiveProps(nextProps) {
      var _state3 = this.state,
          mouseX = _state3.mouseX,
          mouseY = _state3.mouseY,
          resizeFactorX = _state3.resizeFactorX,
          resizeFactorY = _state3.resizeFactorY;
      var _props3 = this.props,
          projection = _props3.projection,
          center = _props3.center,
          zoom = _props3.zoom;


      var zoomFactor = nextProps.zoom / zoom;
      var centerChanged = JSON.stringify(nextProps.center) !== JSON.stringify(center);

      this.setState({
        zoom: nextProps.zoom,
        mouseX: centerChanged ? (0, _utils.calculateMousePosition)("x", nextProps.projection, nextProps, nextProps.zoom, resizeFactorX) : mouseX * zoomFactor,
        mouseY: centerChanged ? (0, _utils.calculateMousePosition)("y", nextProps.projection, nextProps, nextProps.zoom, resizeFactorY) : mouseY * zoomFactor
      });
    }
  }, {
    key: "handleResize",
    value: function handleResize() {
      var _props4 = this.props,
          width = _props4.width,
          height = _props4.height,
          projection = _props4.projection,
          zoom = _props4.zoom;


      var resizeFactorX = (0, _utils.calculateResizeFactor)(this.zoomableGroupNode.parentNode.getBoundingClientRect().width, width);
      var resizeFactorY = (0, _utils.calculateResizeFactor)(this.zoomableGroupNode.parentNode.getBoundingClientRect().height, height);

      var xPercentageChange = 1 / resizeFactorX * this.state.resizeFactorX;
      var yPercentageChange = 1 / resizeFactorY * this.state.resizeFactorY;

      this.setState({
        resizeFactorX: resizeFactorX,
        resizeFactorY: resizeFactorY,
        mouseX: this.state.mouseX * xPercentageChange,
        mouseY: this.state.mouseY * yPercentageChange
      });
    }
  }, {
    key: "componentDidMount",
    value: function componentDidMount() {
      var _props5 = this.props,
          width = _props5.width,
          height = _props5.height,
          projection = _props5.projection,
          zoom = _props5.zoom;


      var resizeFactorX = (0, _utils.calculateResizeFactor)(this.zoomableGroupNode.parentNode.getBoundingClientRect().width, width);
      var resizeFactorY = (0, _utils.calculateResizeFactor)(this.zoomableGroupNode.parentNode.getBoundingClientRect().height, height);

      this.setState({
        resizeFactorX: resizeFactorX,
        resizeFactorY: resizeFactorY,
        mouseX: (0, _utils.calculateMousePosition)("x", projection, this.props, zoom, resizeFactorX),
        mouseY: (0, _utils.calculateMousePosition)("y", projection, this.props, zoom, resizeFactorY)
      });

      window.addEventListener("resize", this.handleResize);
      window.addEventListener("mouseup", this.handleMouseUp);
      this.zoomableGroupNode.addEventListener("touchmove", this.preventTouchScroll);
    }
  }, {
    key: "componentWillUnmount",
    value: function componentWillUnmount() {
      window.removeEventListener("resize", this.handleResize);
      window.removeEventListener("mouseup", this.handleMouseUp);
      this.zoomableGroupNode.removeEventListener("touchmove", this.preventTouchScroll);
    }
  }, {
    key: "render",
    value: function render() {
      var _this2 = this;

      var _props6 = this.props,
          width = _props6.width,
          height = _props6.height,
          zoom = _props6.zoom,
          style = _props6.style,
          projection = _props6.projection,
          children = _props6.children;
      var _state4 = this.state,
          mouseX = _state4.mouseX,
          mouseY = _state4.mouseY,
          resizeFactorX = _state4.resizeFactorX,
          resizeFactorY = _state4.resizeFactorY;


      return _react2.default.createElement(
        "g",
        { className: "rsm-zoomable-group",
          ref: function ref(zoomableGroupNode) {
            return _this2.zoomableGroupNode = zoomableGroupNode;
          },
          transform: "\n           translate(\n             " + Math.round((width / 2 + resizeFactorX * mouseX) * 100) / 100 + "\n             " + Math.round((height / 2 + resizeFactorY * mouseY) * 100) / 100 + "\n           )\n           scale(" + zoom + ")\n           translate(" + -width / 2 + " " + -height / 2 + ")\n         ",
          onMouseMove: this.handleMouseMove,
          onMouseUp: this.handleMouseUp,
          onMouseDown: this.handleMouseDown,
          onTouchStart: this.handleTouchStart,
          onTouchMove: this.handleTouchMove,
          onTouchEnd: this.handleMouseUp,
          style: style
        },
        _react2.default.createElement("rect", {
          x: this.state.backdrop.x,
          y: this.state.backdrop.y,
          width: this.state.backdrop.width,
          height: this.state.backdrop.height,
          fill: "transparent",
          style: { strokeWidth: 0 }
        }),
        (0, _utils.createNewChildren)(children, this.props)
      );
    }
  }]);

  return ZoomableGroup;
}(_react.Component);

ZoomableGroup.defaultProps = {
  center: [0, 0],
  backdrop: {
    x: [-179.9, 179.9],
    y: [89.9, -89.9]
  },
  zoom: 1,
  disablePanning: false
};

exports.default = ZoomableGroup;

/***/ }),
/* 11 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

var _utils = __webpack_require__(5);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var ZoomableGlobe = function (_Component) {
  _inherits(ZoomableGlobe, _Component);

  function ZoomableGlobe(props) {
    _classCallCheck(this, ZoomableGlobe);

    var _this = _possibleConstructorReturn(this, (ZoomableGlobe.__proto__ || Object.getPrototypeOf(ZoomableGlobe)).call(this, props));

    var initialRotation = props.projection.rotate();

    _this.state = {
      mouseX: 0,
      mouseY: 0,
      mouseXStart: 0,
      mouseYStart: 0,
      isPressed: false,
      rotation: [initialRotation[0] - props.center[0], initialRotation[1] - props.center[1], initialRotation[2]]
    };

    _this.handleMouseMove = _this.handleMouseMove.bind(_this);
    _this.handleMouseUp = _this.handleMouseUp.bind(_this);
    _this.handleMouseDown = _this.handleMouseDown.bind(_this);
    _this.handleTouchStart = _this.handleTouchStart.bind(_this);
    _this.handleTouchMove = _this.handleTouchMove.bind(_this);
    return _this;
  }

  _createClass(ZoomableGlobe, [{
    key: "handleMouseMove",
    value: function handleMouseMove(_ref) {
      var pageX = _ref.pageX,
          pageY = _ref.pageY,
          clientX = _ref.clientX,
          clientY = _ref.clientY;

      if (this.props.disablePanning) return;
      if (!this.state.isPressed) return;

      var differenceX = clientX - this.state.mouseXStart;
      var differenceY = clientY - this.state.mouseYStart;

      this.setState({
        mouseX: clientX,
        mouseY: clientY,
        mouseXStart: clientX,
        mouseYStart: clientY,
        rotation: [this.state.rotation[0] + differenceX * this.props.sensitivity, this.state.rotation[1] - differenceY * this.props.sensitivity, this.state.rotation[2]]
      });
    }
  }, {
    key: "handleTouchMove",
    value: function handleTouchMove(_ref2) {
      var touches = _ref2.touches;

      this.handleMouseMove(touches[0]);
    }
  }, {
    key: "handleMouseUp",
    value: function handleMouseUp() {
      if (this.props.disablePanning) return;
      if (!this.state.isPressed) return;
      this.setState({
        isPressed: false
      });
      if (!this.props.onMoveEnd) return;
      var newCenter = this.props.projection.invert([this.props.width / 2, this.props.height / 2]);
      this.props.onMoveEnd(newCenter);
    }
  }, {
    key: "handleMouseDown",
    value: function handleMouseDown(_ref3) {
      var pageX = _ref3.pageX,
          pageY = _ref3.pageY,
          clientX = _ref3.clientX,
          clientY = _ref3.clientY;

      if (this.props.disablePanning) return;
      this.setState({
        isPressed: true,
        mouseXStart: clientX,
        mouseYStart: clientY
      });
      if (!this.props.onMoveStart) return;
      var currentCenter = this.props.projection.invert([this.props.width / 2, this.props.height / 2]);
      this.props.onMoveStart(currentCenter);
    }
  }, {
    key: "handleTouchStart",
    value: function handleTouchStart(_ref4) {
      var touches = _ref4.touches;

      if (touches.length > 1) {
        this.handleMouseDown(touches[0]);
      } else {
        this.handleMouseUp();
      }
    }
  }, {
    key: "preventTouchScroll",
    value: function preventTouchScroll(evt) {
      if (evt.touches.length > 1) {
        evt.preventDefault();
      }
    }
  }, {
    key: "componentWillReceiveProps",
    value: function componentWillReceiveProps(nextProps) {
      var _state = this.state,
          mouseX = _state.mouseX,
          mouseY = _state.mouseY;
      var _props = this.props,
          projection = _props.projection,
          center = _props.center,
          zoom = _props.zoom;


      var zoomFactor = nextProps.zoom / zoom;
      var centerChanged = JSON.stringify(nextProps.center) !== JSON.stringify(center);

      this.setState({
        zoom: nextProps.zoom,
        rotation: centerChanged ? [-nextProps.center[0], -nextProps.center[1], this.state.rotation[2]] : this.state.rotation
      });
    }
  }, {
    key: "componentDidMount",
    value: function componentDidMount() {
      var _props2 = this.props,
          width = _props2.width,
          height = _props2.height,
          projection = _props2.projection,
          zoom = _props2.zoom;


      window.addEventListener("resize", this.handleResize);
      window.addEventListener("mouseup", this.handleMouseUp);
      this.zoomableGlobeNode.addEventListener("touchmove", this.preventTouchScroll);
    }
  }, {
    key: "componentWillUnmount",
    value: function componentWillUnmount() {
      window.removeEventListener("resize", this.handleResize);
      window.removeEventListener("mouseup", this.handleMouseUp);
      this.zoomableGlobeNode.removeEventListener("touchmove", this.preventTouchScroll);
    }
  }, {
    key: "render",
    value: function render() {
      var _this2 = this;

      var _props3 = this.props,
          width = _props3.width,
          height = _props3.height,
          zoom = _props3.zoom,
          style = _props3.style,
          projection = _props3.projection,
          children = _props3.children;
      var _state2 = this.state,
          mouseX = _state2.mouseX,
          mouseY = _state2.mouseY;


      return _react2.default.createElement(
        "g",
        { className: "rsm-zoomable-globe",
          ref: function ref(zoomableGlobeNode) {
            return _this2.zoomableGlobeNode = zoomableGlobeNode;
          },
          transform: "\n           translate(" + width / 2 + " " + height / 2 + ")\n           scale(" + zoom + ")\n           translate(" + -width / 2 + " " + -height / 2 + ")\n         ",
          onMouseMove: this.handleMouseMove,
          onMouseUp: this.handleMouseUp,
          onMouseDown: this.handleMouseDown,
          onTouchStart: this.handleTouchStart,
          onTouchMove: this.handleTouchMove,
          onTouchEnd: this.handleMouseUp,
          style: style
        },
        (0, _utils.createNewChildren)(children, {
          width: width,
          height: height,
          center: this.center,
          backdrop: this.backdrop,
          zoom: this.props.zoom,
          disablePanning: this.props.disablePanning,
          children: children,
          projection: projection.rotate(this.state.rotation)
        })
      );
    }
  }]);

  return ZoomableGlobe;
}(_react.Component);

ZoomableGlobe.defaultProps = {
  center: [0, 0],
  zoom: 1,
  disablePanning: false,
  sensitivity: 0.25
};

exports.default = ZoomableGlobe;

/***/ }),
/* 12 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _topojsonClient = __webpack_require__(22);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var Geographies = function (_Component) {
  _inherits(Geographies, _Component);

  function Geographies(props) {
    _classCallCheck(this, Geographies);

    var _this = _possibleConstructorReturn(this, (Geographies.__proto__ || Object.getPrototypeOf(Geographies)).call(this, props));

    _this.state = {
      geographyPaths: _this.shouldFetchGeographies(props.geography) ? [] : _this.parseGeographies(props.geography)
    };
    return _this;
  }

  _createClass(Geographies, [{
    key: "componentWillReceiveProps",
    value: function componentWillReceiveProps(nextProps) {
      if (nextProps.geography !== this.props.geography) {
        if (this.shouldFetchGeographies(nextProps.geography)) {
          this.fetchGeographies(nextProps.geography);
        } else {
          this.setState({
            geographyPaths: this.parseGeographies(nextProps.geography)
          });
        }
      }
    }
  }, {
    key: "shouldComponentUpdate",
    value: function shouldComponentUpdate(nextProps, nextState) {
      var geoPathsChanged = nextState.geographyPaths.length !== this.state.geographyPaths.length;
      return geoPathsChanged || nextProps.disableOptimization;
    }
  }, {
    key: "componentDidMount",
    value: function componentDidMount() {
      if (this.props.geographyUrl || this.props.geographyPaths) {
        console.warn("You are using the deprecated geographyUrl or geographyPaths props. Use the new geography prop instead. Check out the new docs here: https://github.com/zcreativelabs/react-simple-maps#Geographies-component");
      }
      if (this.shouldFetchGeographies(this.props.geography)) {
        this.fetchGeographies(this.props.geography);
      }
    }
  }, {
    key: "componentWillUnmount",
    value: function componentWillUnmount() {
      this.cancelPendingRequest();
    }
  }, {
    key: "render",
    value: function render() {
      var _props = this.props,
          projection = _props.projection,
          style = _props.style,
          children = _props.children;

      return _react2.default.createElement(
        "g",
        { className: "rsm-geographies", style: style },
        children(this.state.geographyPaths || [], projection)
      );
    }
  }, {
    key: "shouldFetchGeographies",
    value: function shouldFetchGeographies(geography) {
      return typeof geography === 'string';
    }
  }, {
    key: "parseGeographies",
    value: function parseGeographies(geography) {
      if (Array.isArray(geography)) {
        return geography;
      }

      if (Object.prototype.toString.call(geography) === '[object Object]') {
        return (0, _topojsonClient.feature)(geography, geography.objects[Object.keys(geography.objects)[0]]).features;
      }

      return [];
    }
  }, {
    key: "fetchGeographies",
    value: function fetchGeographies(geography) {
      var _this2 = this;

      var request = new XMLHttpRequest();
      request.open("GET", geography, true);
      request.onload = function () {
        if (request.status >= 200 && request.status < 400) {
          var geographyPaths = JSON.parse(request.responseText);
          _this2.setState({
            geographyPaths: _this2.parseGeographies(geographyPaths)
          }, function () {
            if (_this2.props.onGeographyPathsLoaded) {
              _this2.props.onGeographyPathsLoaded(String(request.status));
            }
          });
        } else {
          if (_this2.props.onGeographyPathsLoaded) {
            _this2.props.onGeographyPathsLoaded(String(request.status));
          }
        }
      };
      request.onerror = function () {
        console.log("There was a connection error...");
      };
      request.send();

      this.cancelPendingRequest();
      this._xhr = request;
    }
  }, {
    key: "cancelPendingRequest",
    value: function cancelPendingRequest() {
      if (this._xhr) {
        this._xhr.abort();
        this._xhr = null;
      }
    }
  }]);

  return Geographies;
}(_react.Component);

Geographies.defaultProps = {
  componentIdentifier: "Geographies",
  disableOptimization: false,
  geography: ""
};

exports.default = Geographies;

/***/ }),
/* 13 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; };

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

var _utils = __webpack_require__(5);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var pathCache = {};

var renderPath = function renderPath(cacheId, geography, projection, round, precision) {
  if (pathCache[cacheId]) return pathCache[cacheId];

  var pathString = cacheId ? pathCache[cacheId] ? pathCache[cacheId] : round ? (0, _utils.roundPath)((0, _d3Geo.geoPath)().projection(projection)(geography), precision) : (0, _d3Geo.geoPath)().projection(projection)(geography) : round ? (0, _utils.roundPath)((0, _d3Geo.geoPath)().projection(projection)(geography), precision) : (0, _d3Geo.geoPath)().projection(projection)(geography);

  if (cacheId) pathCache[cacheId] = pathString;

  return pathString;
};

var Geography = function (_Component) {
  _inherits(Geography, _Component);

  function Geography() {
    _classCallCheck(this, Geography);

    var _this = _possibleConstructorReturn(this, (Geography.__proto__ || Object.getPrototypeOf(Geography)).call(this));

    _this.state = {
      hover: false,
      pressed: false
    };

    _this.handleMouseEnter = _this.handleMouseEnter.bind(_this);
    _this.handleMouseMove = _this.handleMouseMove.bind(_this);
    _this.handleMouseLeave = _this.handleMouseLeave.bind(_this);
    _this.handleMouseDown = _this.handleMouseDown.bind(_this);
    _this.handleMouseUp = _this.handleMouseUp.bind(_this);
    _this.handleMouseClick = _this.handleMouseClick.bind(_this);
    _this.handleFocus = _this.handleFocus.bind(_this);
    _this.handleBlur = _this.handleBlur.bind(_this);
    return _this;
  }

  _createClass(Geography, [{
    key: "handleMouseClick",
    value: function handleMouseClick(evt) {
      evt.persist();
      var _props = this.props,
          onClick = _props.onClick,
          geography = _props.geography;

      return onClick && onClick(geography, evt);
    }
  }, {
    key: "handleMouseEnter",
    value: function handleMouseEnter(evt) {
      evt.persist();
      var _props2 = this.props,
          onMouseEnter = _props2.onMouseEnter,
          geography = _props2.geography;

      this.setState({
        hover: true
      }, function () {
        return onMouseEnter && onMouseEnter(geography, evt);
      });
    }
  }, {
    key: "handleMouseMove",
    value: function handleMouseMove(evt) {
      evt.persist();
      if (this.state.pressed) return;
      var _props3 = this.props,
          onMouseMove = _props3.onMouseMove,
          geography = _props3.geography;

      if (!this.state.hover) {
        this.setState({
          hover: true
        }, function () {
          return onMouseMove && onMouseMove(geography, evt);
        });
      } else if (onMouseMove) onMouseMove(geography, evt);else return;
    }
  }, {
    key: "handleMouseLeave",
    value: function handleMouseLeave(evt) {
      evt.persist();
      var _props4 = this.props,
          onMouseLeave = _props4.onMouseLeave,
          geography = _props4.geography;

      this.setState({
        hover: false,
        pressed: false
      }, function () {
        return onMouseLeave && onMouseLeave(geography, evt);
      });
    }
  }, {
    key: "handleMouseDown",
    value: function handleMouseDown(evt) {
      evt.persist();
      var _props5 = this.props,
          onMouseDown = _props5.onMouseDown,
          geography = _props5.geography;

      this.setState({
        pressed: true
      }, function () {
        return onMouseDown && onMouseDown(geography, evt);
      });
    }
  }, {
    key: "handleMouseUp",
    value: function handleMouseUp(evt) {
      evt.persist();
      var _props6 = this.props,
          onMouseUp = _props6.onMouseUp,
          geography = _props6.geography;

      this.setState({
        pressed: false
      }, function () {
        return onMouseUp && onMouseUp(geography, evt);
      });
    }
  }, {
    key: "handleFocus",
    value: function handleFocus(evt) {
      evt.persist();
      var _props7 = this.props,
          onFocus = _props7.onFocus,
          geography = _props7.geography;

      this.setState({
        hover: true
      }, function () {
        return onFocus && onFocus(geography, evt);
      });
    }
  }, {
    key: "handleBlur",
    value: function handleBlur(evt) {
      evt.persist();
      var _props8 = this.props,
          onBlur = _props8.onBlur,
          geography = _props8.geography;

      this.setState({
        hover: false
      }, function () {
        return onBlur && onBlur(geography, evt);
      });
    }
  }, {
    key: "render",
    value: function render() {
      var _this2 = this;

      var _props9 = this.props,
          geography = _props9.geography,
          projection = _props9.projection,
          round = _props9.round,
          cacheId = _props9.cacheId,
          precision = _props9.precision,
          tabable = _props9.tabable,
          style = _props9.style;
      var _state = this.state,
          hover = _state.hover,
          pressed = _state.pressed;


      var pathString = renderPath(cacheId, geography, projection, round, precision);

      var excludeProps = ["geography", "projection", "round", "cacheId", "precision", "tabable", "style", "onClick", "onMouseEnter", "onMouseMove", "onMouseLeave", "onMouseDown", "onMouseUp", "onFocus", "onBlur"];

      var restProps = Object.keys(this.props).filter(function (key) {
        return excludeProps.indexOf(key) === -1;
      }).reduce(function (obj, key) {
        obj[key] = _this2.props[key];
        return obj;
      }, {});

      return _react2.default.createElement("path", _extends({
        d: pathString,
        className: "rsm-geography" + (pressed ? " rsm-geography--pressed" : "") + (hover ? " rsm-geography--hover" : ""),
        style: style[pressed || hover ? pressed ? "pressed" : "hover" : "default"],
        onClick: this.handleMouseClick,
        onMouseEnter: this.handleMouseEnter,
        onMouseMove: this.handleMouseMove,
        onMouseLeave: this.handleMouseLeave,
        onMouseDown: this.handleMouseDown,
        onMouseUp: this.handleMouseUp,
        onFocus: this.handleFocus,
        onBlur: this.handleBlur,
        tabIndex: tabable ? 0 : -1
      }, restProps));
    }
  }]);

  return Geography;
}(_react.Component);

Geography.defaultProps = {
  precision: 0.1,
  cacheId: null,
  round: false,
  tabable: true,
  style: {
    default: {},
    hover: {},
    pressed: {}
  }
};

exports.default = Geography;

/***/ }),
/* 14 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var Marker = function (_Component) {
  _inherits(Marker, _Component);

  function Marker() {
    _classCallCheck(this, Marker);

    var _this = _possibleConstructorReturn(this, (Marker.__proto__ || Object.getPrototypeOf(Marker)).call(this));

    _this.state = {
      hover: false,
      pressed: false
    };

    _this.handleMouseEnter = _this.handleMouseEnter.bind(_this);
    _this.handleMouseLeave = _this.handleMouseLeave.bind(_this);
    _this.handleMouseDown = _this.handleMouseDown.bind(_this);
    _this.handleMouseUp = _this.handleMouseUp.bind(_this);
    _this.handleMouseClick = _this.handleMouseClick.bind(_this);
    _this.handleMouseMove = _this.handleMouseMove.bind(_this);
    _this.handleFocus = _this.handleFocus.bind(_this);
    _this.handleBlur = _this.handleBlur.bind(_this);
    return _this;
  }

  _createClass(Marker, [{
    key: "handleMouseEnter",
    value: function handleMouseEnter(evt) {
      evt.persist();
      var _props = this.props,
          onMouseEnter = _props.onMouseEnter,
          marker = _props.marker;

      this.setState({
        hover: true
      }, function () {
        return onMouseEnter && onMouseEnter(marker, evt);
      });
    }
  }, {
    key: "handleMouseMove",
    value: function handleMouseMove(evt) {
      evt.persist();
      if (this.state.pressed) return;
      var _props2 = this.props,
          onMouseMove = _props2.onMouseMove,
          marker = _props2.marker;

      if (!this.state.hover) {
        this.setState({
          hover: true
        }, function () {
          return onMouseMove && onMouseMove(marker, evt);
        });
      } else if (onMouseMove) onMouseMove(marker, evt);else return;
    }
  }, {
    key: "handleMouseLeave",
    value: function handleMouseLeave(evt) {
      evt.persist();
      var _props3 = this.props,
          onMouseLeave = _props3.onMouseLeave,
          marker = _props3.marker;

      this.setState({
        hover: false
      }, function () {
        return onMouseLeave && onMouseLeave(marker, evt);
      });
    }
  }, {
    key: "handleMouseDown",
    value: function handleMouseDown(evt) {
      evt.persist();
      var _props4 = this.props,
          onMouseDown = _props4.onMouseDown,
          marker = _props4.marker;

      this.setState({
        pressed: true
      }, function () {
        return onMouseDown && onMouseDown(marker, evt);
      });
    }
  }, {
    key: "handleMouseUp",
    value: function handleMouseUp(evt) {
      evt.persist();
      var _props5 = this.props,
          onMouseUp = _props5.onMouseUp,
          marker = _props5.marker;

      this.setState({
        pressed: false
      }, function () {
        return onMouseUp && onMouseUp(marker, evt);
      });
    }
  }, {
    key: "handleMouseClick",
    value: function handleMouseClick(evt) {
      if (!this.props.onClick) return;
      evt.persist();
      var _props6 = this.props,
          onClick = _props6.onClick,
          marker = _props6.marker,
          projection = _props6.projection;

      return onClick && onClick(marker, projection(marker.coordinates), evt);
    }
  }, {
    key: "handleFocus",
    value: function handleFocus(evt) {
      evt.persist();
      var _props7 = this.props,
          onFocus = _props7.onFocus,
          marker = _props7.marker;

      this.setState({
        hover: true
      }, function () {
        return onFocus && onFocus(marker, evt);
      });
    }
  }, {
    key: "handleBlur",
    value: function handleBlur(evt) {
      evt.persist();
      var _props8 = this.props,
          onBlur = _props8.onBlur,
          marker = _props8.marker;

      this.setState({
        hover: false
      }, function () {
        return onBlur && onBlur(marker, evt);
      });
    }
  }, {
    key: "render",
    value: function render() {
      var _props9 = this.props,
          projection = _props9.projection,
          marker = _props9.marker,
          style = _props9.style,
          tabable = _props9.tabable,
          zoom = _props9.zoom,
          children = _props9.children,
          preserveMarkerAspect = _props9.preserveMarkerAspect,
          width = _props9.width,
          height = _props9.height;
      var _state = this.state,
          pressed = _state.pressed,
          hover = _state.hover;


      var scale = preserveMarkerAspect ? " scale(" + 1 / zoom + ")" : "";
      var translation = projection(marker.coordinates);

      var lineString = {
        "type": "Feature",
        "geometry": {
          "type": "LineString",
          "coordinates": [projection.invert([width / 2, height / 2]), marker.coordinates]
        }
      };

      var radians = Math.PI / 2,
          degrees = 90;
      var isGlobe = projection.clipAngle && projection.clipAngle() === degrees;
      var isHidden = isGlobe && (0, _d3Geo.geoLength)(lineString) > radians;

      return _react2.default.createElement(
        "g",
        { className: "rsm-marker" + (pressed ? " rsm-marker--pressed" : "") + (hover ? " rsm-marker--hover" : ""),
          transform: "translate(\n           " + translation[0] + "\n           " + translation[1] + "\n         ) " + scale,
          style: style[isHidden ? "hidden" : pressed || hover ? pressed ? "pressed" : "hover" : "default"],
          onMouseEnter: this.handleMouseEnter,
          onMouseLeave: this.handleMouseLeave,
          onMouseDown: this.handleMouseDown,
          onMouseUp: this.handleMouseUp,
          onClick: this.handleMouseClick,
          onMouseMove: this.handleMouseMove,
          onFocus: this.handleFocus,
          onBlur: this.handleBlur,
          tabIndex: tabable ? 0 : -1
        },
        children
      );
    }
  }]);

  return Marker;
}(_react.Component);

Marker.defaultProps = {
  style: {
    default: {},
    hover: {},
    pressed: {}
  },
  marker: {
    coordinates: [0, 0]
  },
  tabable: true,
  preserveMarkerAspect: true
};

exports.default = Marker;

/***/ }),
/* 15 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; };

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _MapGroup = __webpack_require__(6);

var _MapGroup2 = _interopRequireDefault(_MapGroup);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _objectWithoutProperties(obj, keys) { var target = {}; for (var i in obj) { if (keys.indexOf(i) >= 0) continue; if (!Object.prototype.hasOwnProperty.call(obj, i)) continue; target[i] = obj[i]; } return target; }

var Markers = function Markers(_ref) {
  var groupName = _ref.groupName,
      itemName = _ref.itemName,
      componentIdentifier = _ref.componentIdentifier,
      restProps = _objectWithoutProperties(_ref, ["groupName", "itemName", "componentIdentifier"]);

  return _react2.default.createElement(_MapGroup2.default, _extends({
    groupName: groupName,
    itemName: itemName
  }, restProps));
};

Markers.defaultProps = {
  componentIdentifier: "Markers",
  groupName: "markers",
  itemName: "marker"
};

exports.default = Markers;

/***/ }),
/* 16 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var Line = function (_Component) {
  _inherits(Line, _Component);

  function Line(props) {
    _classCallCheck(this, Line);

    var _this = _possibleConstructorReturn(this, (Line.__proto__ || Object.getPrototypeOf(Line)).call(this, props));

    _this.state = {
      hover: false,
      pressed: false
    };

    _this.handleMouseEnter = _this.handleMouseEnter.bind(_this);
    _this.handleMouseLeave = _this.handleMouseLeave.bind(_this);
    _this.handleMouseDown = _this.handleMouseDown.bind(_this);
    _this.handleMouseUp = _this.handleMouseUp.bind(_this);
    _this.handleMouseClick = _this.handleMouseClick.bind(_this);
    _this.handleMouseMove = _this.handleMouseMove.bind(_this);
    _this.handleFocus = _this.handleFocus.bind(_this);
    _this.handleBlur = _this.handleBlur.bind(_this);
    return _this;
  }

  _createClass(Line, [{
    key: "handleMouseEnter",
    value: function handleMouseEnter(evt) {
      evt.persist();
      var _props = this.props,
          onMouseEnter = _props.onMouseEnter,
          line = _props.line;

      this.setState({
        hover: true
      }, function () {
        return onMouseEnter && onMouseEnter(line, evt);
      });
    }
  }, {
    key: "handleMouseMove",
    value: function handleMouseMove(evt) {
      evt.persist();
      if (this.state.pressed) return;
      var _props2 = this.props,
          onMouseMove = _props2.onMouseMove,
          line = _props2.line;

      if (!this.state.hover) {
        this.setState({
          hover: true
        }, function () {
          return onMouseMove && onMouseMove(line, evt);
        });
      } else if (onMouseMove) onMouseMove(line, evt);else return;
    }
  }, {
    key: "handleMouseLeave",
    value: function handleMouseLeave(evt) {
      evt.persist();
      var _props3 = this.props,
          onMouseLeave = _props3.onMouseLeave,
          line = _props3.line;

      this.setState({
        hover: false
      }, function () {
        return onMouseLeave && onMouseLeave(line, evt);
      });
    }
  }, {
    key: "handleMouseDown",
    value: function handleMouseDown(evt) {
      evt.persist();
      var _props4 = this.props,
          onMouseDown = _props4.onMouseDown,
          line = _props4.line;

      this.setState({
        pressed: true
      }, function () {
        return onMouseDown && onMouseDown(line, evt);
      });
    }
  }, {
    key: "handleMouseUp",
    value: function handleMouseUp(evt) {
      evt.persist();
      var _props5 = this.props,
          onMouseUp = _props5.onMouseUp,
          line = _props5.line;

      this.setState({
        pressed: false
      }, function () {
        return onMouseUp && onMouseUp(line, evt);
      });
    }
  }, {
    key: "handleMouseClick",
    value: function handleMouseClick(evt) {
      if (!this.props.onClick) return;
      evt.persist();
      var _props6 = this.props,
          onClick = _props6.onClick,
          line = _props6.line,
          projection = _props6.projection;

      return onClick && onClick(line, [projection(line.coordinates.start), projection(line.coordinates.end)], evt);
    }
  }, {
    key: "handleFocus",
    value: function handleFocus(evt) {
      evt.persist();
      var _props7 = this.props,
          onFocus = _props7.onFocus,
          line = _props7.line;

      this.setState({
        hover: true
      }, function () {
        return onFocus && onFocus(line, evt);
      });
    }
  }, {
    key: "handleBlur",
    value: function handleBlur(evt) {
      evt.persist();
      var _props8 = this.props,
          onBlur = _props8.onBlur,
          line = _props8.line;

      this.setState({
        hover: false
      }, function () {
        return onBlur && onBlur(line, evt);
      });
    }
  }, {
    key: "render",
    value: function render() {
      var _props9 = this.props,
          className = _props9.className,
          projection = _props9.projection,
          line = _props9.line,
          style = _props9.style,
          tabable = _props9.tabable,
          zoom = _props9.zoom,
          preserveMarkerAspect = _props9.preserveMarkerAspect,
          width = _props9.width,
          height = _props9.height,
          buildPath = _props9.buildPath,
          strokeWidth = _props9.strokeWidth;
      var _state = this.state,
          pressed = _state.pressed,
          hover = _state.hover;


      var scale = preserveMarkerAspect ? " scale(" + 1 / zoom + ")" : "";

      var buildLineString = function buildLineString(coordinates) {
        return {
          type: "Feature",
          geometry: {
            type: "LineString",
            coordinates: [projection.invert([width / 2, height / 2]), coordinates]
          }
        };
      };
      var startLineString = buildLineString(line.coordinates.start);
      var endLineString = buildLineString(line.coordinates.end);

      var radians = Math.PI / 2,
          degrees = 90;
      var isGlobe = projection.clipAngle && projection.clipAngle() === degrees;
      var isHidden = isGlobe && ((0, _d3Geo.geoLength)(startLineString) > radians || (0, _d3Geo.geoLength)(endLineString) > radians);

      var start = projection(line.coordinates.start);
      var end = projection(line.coordinates.end);

      var path = buildPath ? buildPath(start, end, line) : "M " + start.join(" ") + " L " + end.join(" ");

      return _react2.default.createElement("path", {
        className: "rsm-line" + (pressed ? " rsm-line--pressed" : "") + (hover ? " rsm-line--hover" : "") + " " + className,
        transform: "" + scale,
        style: style[isHidden ? "hidden" : pressed || hover ? pressed ? "pressed" : "hover" : "default"],
        onMouseEnter: this.handleMouseEnter,
        onMouseLeave: this.handleMouseLeave,
        onMouseDown: this.handleMouseDown,
        onMouseUp: this.handleMouseUp,
        onClick: this.handleMouseClick,
        onMouseMove: this.handleMouseMove,
        onFocus: this.handleFocus,
        onBlur: this.handleBlur,
        tabIndex: tabable ? 0 : -1,
        d: path,
        strokeWidth: strokeWidth
      });
    }
  }]);

  return Line;
}(_react.Component);

Line.defaultProps = {
  style: {
    default: {},
    hover: {},
    pressed: {}
  },
  line: {
    coordinates: {
      start: [0, 0],
      end: [-99.1, 19.4]
    }
  },
  tabable: true,
  preserveMarkerAspect: true,
  strokeWidth: 1,
  className: ""
};

exports.default = Line;

/***/ }),
/* 17 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; };

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _MapGroup = __webpack_require__(6);

var _MapGroup2 = _interopRequireDefault(_MapGroup);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _objectWithoutProperties(obj, keys) { var target = {}; for (var i in obj) { if (keys.indexOf(i) >= 0) continue; if (!Object.prototype.hasOwnProperty.call(obj, i)) continue; target[i] = obj[i]; } return target; }

var Lines = function Lines(_ref) {
  var groupName = _ref.groupName,
      itemName = _ref.itemName,
      componentIdentifier = _ref.componentIdentifier,
      restProps = _objectWithoutProperties(_ref, ["groupName", "itemName", "componentIdentifier"]);

  return _react2.default.createElement(_MapGroup2.default, _extends({
    groupName: groupName,
    itemName: itemName
  }, restProps));
};

Lines.defaultProps = {
  componentIdentifier: "Lines",
  groupName: "lines",
  itemName: "line"
};

exports.default = Lines;

/***/ }),
/* 18 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; };

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

var _utils = __webpack_require__(5);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var Annotation = function (_Component) {
  _inherits(Annotation, _Component);

  function Annotation() {
    _classCallCheck(this, Annotation);

    return _possibleConstructorReturn(this, (Annotation.__proto__ || Object.getPrototypeOf(Annotation)).apply(this, arguments));
  }

  _createClass(Annotation, [{
    key: "render",
    value: function render() {
      var _props = this.props,
          projection = _props.projection,
          subject = _props.subject,
          style = _props.style,
          hiddenStyle = _props.hiddenStyle,
          dx = _props.dx,
          dy = _props.dy,
          zoom = _props.zoom,
          stroke = _props.stroke,
          strokeWidth = _props.strokeWidth,
          children = _props.children,
          curve = _props.curve,
          markerEnd = _props.markerEnd,
          width = _props.width,
          height = _props.height;


      var connectorPath = (0, _utils.createConnectorPath)(null, [-dx / zoom, -dy / zoom], curve);
      var translation = projection(subject);

      var lineString = {
        "type": "Feature",
        "geometry": {
          "type": "LineString",
          "coordinates": [projection.invert([width / 2, height / 2]), subject]
        }
      };

      var radians = Math.PI / 2,
          degrees = 90;
      var isGlobe = projection.clipAngle && projection.clipAngle() === degrees;
      var isHidden = isGlobe && (0, _d3Geo.geoLength)(lineString) > radians;

      return _react2.default.createElement(
        "g",
        {
          className: "rsm-annotation",
          style: isHidden ? _extends({}, style, hiddenStyle) : style,
          transform: "translate(\n          " + (translation[0] + dx / zoom) + "\n          " + (translation[1] + dy / zoom) + "\n        )",
          textAnchor: (0, _utils.createTextAnchor)(dx)
        },
        children,
        _react2.default.createElement("path", {
          d: connectorPath,
          stroke: stroke,
          strokeWidth: strokeWidth,
          fill: "none",
          markerEnd: markerEnd
        })
      );
    }
  }]);

  return Annotation;
}(_react.Component);

Annotation.defaultProps = {
  curve: 0,
  markerEnd: "none",
  componentIdentifier: "Annotation",
  stroke: "#000000",
  strokeWidth: 1,
  zoom: 1
};

exports.default = Annotation;

/***/ }),
/* 19 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; };

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _MapGroup = __webpack_require__(6);

var _MapGroup2 = _interopRequireDefault(_MapGroup);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _objectWithoutProperties(obj, keys) { var target = {}; for (var i in obj) { if (keys.indexOf(i) >= 0) continue; if (!Object.prototype.hasOwnProperty.call(obj, i)) continue; target[i] = obj[i]; } return target; }

var Annotations = function Annotations(_ref) {
  var groupName = _ref.groupName,
      itemName = _ref.itemName,
      componentIdentifier = _ref.componentIdentifier,
      restProps = _objectWithoutProperties(_ref, ["groupName", "itemName", "componentIdentifier"]);

  return _react2.default.createElement(_MapGroup2.default, _extends({
    groupName: groupName,
    itemName: itemName
  }, restProps));
};

Annotations.defaultProps = {
  componentIdentifier: "Annotations",
  groupName: "annotations",
  itemName: "annotation"
};

exports.default = Annotations;

/***/ }),
/* 20 */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var _createClass = function () { function defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } } return function (Constructor, protoProps, staticProps) { if (protoProps) defineProperties(Constructor.prototype, protoProps); if (staticProps) defineProperties(Constructor, staticProps); return Constructor; }; }();

var _react = __webpack_require__(1);

var _react2 = _interopRequireDefault(_react);

var _d3Geo = __webpack_require__(0);

var _utils = __webpack_require__(5);

function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : subClass.__proto__ = superClass; }

var computeGraticule = function computeGraticule(projection, step) {
  return (0, _d3Geo.geoPath)().projection(projection)((0, _d3Geo.geoGraticule)().step(step)());
};

var computeOutline = function computeOutline(projection) {
  return (0, _d3Geo.geoPath)().projection(projection)((0, _d3Geo.geoGraticule)().outline());
};

var Graticule = function (_Component) {
  _inherits(Graticule, _Component);

  function Graticule() {
    _classCallCheck(this, Graticule);

    var _this = _possibleConstructorReturn(this, (Graticule.__proto__ || Object.getPrototypeOf(Graticule)).call(this));

    _this.state = {
      renderGraticule: false,
      graticulePath: "",
      outlinePath: ""
    };
    _this.renderGraticule = _this.renderGraticule.bind(_this);
    return _this;
  }

  _createClass(Graticule, [{
    key: "componentDidMount",
    value: function componentDidMount() {
      this.renderGraticule();
    }
  }, {
    key: "renderGraticule",
    value: function renderGraticule() {
      var _props = this.props,
          step = _props.step,
          projection = _props.projection,
          round = _props.round,
          precision = _props.precision;


      this.setState({
        renderGraticule: true,
        graticulePath: round ? (0, _utils.roundPath)(computeGraticule(projection, step), precision) : computeGraticule(projection, step),
        outlinePath: round ? (0, _utils.roundPath)(computeOutline(projection), precision) : computeOutline(projection)
      });
    }
  }, {
    key: "componentWillReceiveProps",
    value: function componentWillReceiveProps(nextProps) {
      var _props2 = this.props,
          step = _props2.step,
          projection = _props2.projection,
          round = _props2.round,
          precision = _props2.precision,
          globe = _props2.globe;


      if (nextProps.round !== round || nextProps.precision !== precision || globe) {
        this.setState({
          graticulePath: nextProps.round ? (0, _utils.roundPath)(computeGraticule(projection, step), precision) : computeGraticule(projection, step),
          outlinePath: nextProps.round ? (0, _utils.roundPath)(computeOutline(projection), precision) : computeOutline(projection)
        });
      }
    }
  }, {
    key: "shouldComponentUpdate",
    value: function shouldComponentUpdate(nextProps) {
      return nextProps.disableOptimization;
    }
  }, {
    key: "render",
    value: function render() {
      var _props3 = this.props,
          zoom = _props3.zoom,
          style = _props3.style,
          outline = _props3.outline,
          fill = _props3.fill,
          stroke = _props3.stroke;


      return this.state.renderGraticule && _react2.default.createElement(
        "g",
        { className: "rsm-graticule" },
        _react2.default.createElement("path", {
          fill: fill,
          stroke: stroke,
          d: this.state.graticulePath,
          style: style
        }),
        outline && _react2.default.createElement("path", {
          fill: fill,
          stroke: stroke,
          d: this.state.outlinePath,
          style: style
        })
      );
    }
  }]);

  return Graticule;
}(_react.Component);

Graticule.defaultProps = {
  componentIdentifier: "Graticule",
  disableOptimization: true,
  globe: false,
  round: true,
  precision: 0.1,
  step: [10, 10],
  outline: true,
  stroke: "#DDDDDD",
  fill: "transparent",
  style: {
    pointerEvents: "none"
  }
};

exports.default = Graticule;

/***/ }),
/* 21 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
__webpack_require__.r(__webpack_exports__);

// EXTERNAL MODULE: ./node_modules/d3-geo/index.js + 55 modules
var d3_geo = __webpack_require__(0);

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/math.js
var abs = Math.abs;
var atan = Math.atan;
var atan2 = Math.atan2;
var ceil = Math.ceil;
var cos = Math.cos;
var exp = Math.exp;
var floor = Math.floor;
var log = Math.log;
var max = Math.max;
var min = Math.min;
var pow = Math.pow;
var round = Math.round;
var sign = Math.sign || function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; };
var sin = Math.sin;
var tan = Math.tan;

var math_epsilon = 1e-6;
var epsilon2 = 1e-12;
var pi = Math.PI;
var halfPi = pi / 2;
var quarterPi = pi / 4;
var sqrt1_2 = Math.SQRT1_2;
var sqrt2 = sqrt(2);
var sqrtPi = sqrt(pi);
var tau = pi * 2;
var degrees = 180 / pi;
var radians = pi / 180;

function sinci(x) {
  return x ? x / Math.sin(x) : 1;
}

function asin(x) {
  return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
}

function acos(x) {
  return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
}

function sqrt(x) {
  return x > 0 ? Math.sqrt(x) : 0;
}

function tanh(x) {
  x = exp(2 * x);
  return (x - 1) / (x + 1);
}

function sinh(x) {
  return (exp(x) - exp(-x)) / 2;
}

function cosh(x) {
  return (exp(x) + exp(-x)) / 2;
}

function arsinh(x) {
  return log(x + sqrt(x * x + 1));
}

function arcosh(x) {
  return log(x + sqrt(x * x - 1));
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/airy.js



function airyRaw(beta) {
  var tanBeta_2 = tan(beta / 2),
      b = 2 * log(cos(beta / 2)) / (tanBeta_2 * tanBeta_2);

  function forward(x, y) {
    var cosx = cos(x),
        cosy = cos(y),
        siny = sin(y),
        cosz = cosy * cosx,
        k = -((1 - cosz ? log((1 + cosz) / 2) / (1 - cosz) : -0.5) + b / (1 + cosz));
    return [k * cosy * sin(x), k * siny];
  }

  forward.invert = function(x, y) {
    var r = sqrt(x * x + y * y),
        z = -beta / 2,
        i = 50, delta;
    if (!r) return [0, 0];
    do {
      var z_2 = z / 2,
          cosz_2 = cos(z_2),
          sinz_2 = sin(z_2),
          tanz_2 = tan(z_2),
          lnsecz_2 = log(1 / cosz_2);
      z -= delta = (2 / tanz_2 * lnsecz_2 - b * tanz_2 - r) / (-lnsecz_2 / (sinz_2 * sinz_2) + 1 - b / (2 * cosz_2 * cosz_2));
    } while (abs(delta) > math_epsilon && --i > 0);
    var sinz = sin(z);
    return [atan2(x * sinz, r * cos(z)), asin(y * sinz / r)];
  };

  return forward;
}

/* harmony default export */ var airy = (function() {
  var beta = halfPi,
      m = Object(d3_geo["geoProjectionMutator"])(airyRaw),
      p = m(beta);

  p.radius = function(_) {
    return arguments.length ? m(beta = _ * radians) : beta * degrees;
  };

  return p
      .scale(179.976)
      .clipAngle(147);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/aitoff.js



function aitoffRaw(x, y) {
  var cosy = cos(y), sincia = sinci(acos(cosy * cos(x /= 2)));
  return [2 * cosy * sin(x) * sincia, sin(y) * sincia];
}

// Abort if [x, y] is not within an ellipse centered at [0, 0] with
// semi-major axis pi and semi-minor axis pi/2.
aitoffRaw.invert = function(x, y) {
  if (x * x + 4 * y * y > pi * pi + math_epsilon) return;
  var x1 = x, y1 = y, i = 25;
  do {
    var sinx = sin(x1),
        sinx_2 = sin(x1 / 2),
        cosx_2 = cos(x1 / 2),
        siny = sin(y1),
        cosy = cos(y1),
        sin_2y = sin(2 * y1),
        sin2y = siny * siny,
        cos2y = cosy * cosy,
        sin2x_2 = sinx_2 * sinx_2,
        c = 1 - cos2y * cosx_2 * cosx_2,
        e = c ? acos(cosy * cosx_2) * sqrt(f = 1 / c) : f = 0,
        f,
        fx = 2 * e * cosy * sinx_2 - x,
        fy = e * siny - y,
        dxdx = f * (cos2y * sin2x_2 + e * cosy * cosx_2 * sin2y),
        dxdy = f * (0.5 * sinx * sin_2y - e * 2 * siny * sinx_2),
        dydx = f * 0.25 * (sin_2y * sinx_2 - e * siny * cos2y * sinx),
        dydy = f * (sin2y * cosx_2 + e * sin2x_2 * cosy),
        z = dxdy * dydx - dydy * dxdx;
    if (!z) break;
    var dx = (fy * dxdy - fx * dydy) / z,
        dy = (fx * dydx - fy * dxdx) / z;
    x1 -= dx, y1 -= dy;
  } while ((abs(dx) > math_epsilon || abs(dy) > math_epsilon) && --i > 0);
  return [x1, y1];
};

/* harmony default export */ var aitoff = (function() {
  return Object(d3_geo["geoProjection"])(aitoffRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/armadillo.js



function armadilloRaw(phi0) {
  var sinPhi0 = sin(phi0),
      cosPhi0 = cos(phi0),
      sPhi0 = phi0 >= 0 ? 1 : -1,
      tanPhi0 = tan(sPhi0 * phi0),
      k = (1 + sinPhi0 - cosPhi0) / 2;

  function forward(lambda, phi) {
    var cosPhi = cos(phi),
        cosLambda = cos(lambda /= 2);
    return [
      (1 + cosPhi) * sin(lambda),
      (sPhi0 * phi > -atan2(cosLambda, tanPhi0) - 1e-3 ? 0 : -sPhi0 * 10) + k + sin(phi) * cosPhi0 - (1 + cosPhi) * sinPhi0 * cosLambda // TODO D3 core should allow null or [NaN, NaN] to be returned.
    ];
  }

  forward.invert = function(x, y) {
    var lambda = 0,
        phi = 0,
        i = 50;
    do {
      var cosLambda = cos(lambda),
          sinLambda = sin(lambda),
          cosPhi = cos(phi),
          sinPhi = sin(phi),
          A = 1 + cosPhi,
          fx = A * sinLambda - x,
          fy = k + sinPhi * cosPhi0 - A * sinPhi0 * cosLambda - y,
          dxdLambda = A * cosLambda / 2,
          dxdPhi = -sinLambda * sinPhi,
          dydLambda = sinPhi0 * A * sinLambda / 2,
          dydPhi = cosPhi0 * cosPhi + sinPhi0 * cosLambda * sinPhi,
          denominator = dxdPhi * dydLambda - dydPhi * dxdLambda,
          dLambda = (fy * dxdPhi - fx * dydPhi) / denominator / 2,
          dPhi = (fx * dydLambda - fy * dxdLambda) / denominator;
      lambda -= dLambda, phi -= dPhi;
    } while ((abs(dLambda) > math_epsilon || abs(dPhi) > math_epsilon) && --i > 0);
    return sPhi0 * phi > -atan2(cos(lambda), tanPhi0) - 1e-3 ? [lambda * 2, phi] : null;
  };

  return forward;
}

/* harmony default export */ var armadillo = (function() {
  var phi0 = 20 * radians,
      sPhi0 = phi0 >= 0 ? 1 : -1,
      tanPhi0 = tan(sPhi0 * phi0),
      m = Object(d3_geo["geoProjectionMutator"])(armadilloRaw),
      p = m(phi0),
      stream_ = p.stream;

  p.parallel = function(_) {
    if (!arguments.length) return phi0 * degrees;
    tanPhi0 = tan((sPhi0 = (phi0 = _ * radians) >= 0 ? 1 : -1) * phi0);
    return m(phi0);
  };

  p.stream = function(stream) {
    var rotate = p.rotate(),
        rotateStream = stream_(stream),
        sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart(), sphereStream.lineStart();
      for (var lambda = sPhi0 * -180; sPhi0 * lambda < 180; lambda += sPhi0 * 90) sphereStream.point(lambda, sPhi0 * 90);
      while (sPhi0 * (lambda -= phi0) >= -180) { // TODO precision?
        sphereStream.point(lambda, sPhi0 * -atan2(cos(lambda * radians / 2), tanPhi0) * degrees);
      }
      sphereStream.lineEnd(), sphereStream.polygonEnd();
    };
    return rotateStream;
  };

  return p
      .scale(218.695)
      .center([0, 28.0974]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/august.js



function augustRaw(lambda, phi) {
  var tanPhi = tan(phi / 2),
      k = sqrt(1 - tanPhi * tanPhi),
      c = 1 + k * cos(lambda /= 2),
      x = sin(lambda) * k / c,
      y = tanPhi / c,
      x2 = x * x,
      y2 = y * y;
  return [
    4 / 3 * x * (3 + x2 - 3 * y2),
    4 / 3 * y * (3 + 3 * x2 - y2)
  ];
}

augustRaw.invert = function(x, y) {
  x *= 3 / 8, y *= 3 / 8;
  if (!x && abs(y) > 1) return null;
  var x2 = x * x,
      y2 = y * y,
      s = 1 + x2 + y2,
      sin3Eta = sqrt((s - sqrt(s * s - 4 * y * y)) / 2),
      eta = asin(sin3Eta) / 3,
      xi = sin3Eta ? arcosh(abs(y / sin3Eta)) / 3 : arsinh(abs(x)) / 3,
      cosEta = cos(eta),
      coshXi = cosh(xi),
      d = coshXi * coshXi - cosEta * cosEta;
  return [
    sign(x) * 2 * atan2(sinh(xi) * cosEta, 0.25 - d),
    sign(y) * 2 * atan2(coshXi * sin(eta), 0.25 + d)
  ];
};

/* harmony default export */ var august = (function() {
  return Object(d3_geo["geoProjection"])(augustRaw)
      .scale(66.1603);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/baker.js



var sqrt8 = sqrt(8),
    baker_phi0 = log(1 + sqrt2);

function bakerRaw(lambda, phi) {
  var phi0 = abs(phi);
  return phi0 < quarterPi
      ? [lambda, log(tan(quarterPi + phi / 2))]
      : [lambda * cos(phi0) * (2 * sqrt2 - 1 / sin(phi0)), sign(phi) * (2 * sqrt2 * (phi0 - quarterPi) - log(tan(phi0 / 2)))];
}

bakerRaw.invert = function(x, y) {
  if ((y0 = abs(y)) < baker_phi0) return [x, 2 * atan(exp(y)) - halfPi];
  var phi = quarterPi, i = 25, delta, y0;
  do {
    var cosPhi_2 = cos(phi / 2), tanPhi_2 = tan(phi / 2);
    phi -= delta = (sqrt8 * (phi - quarterPi) - log(tanPhi_2) - y0) / (sqrt8 - cosPhi_2 * cosPhi_2 / (2 * tanPhi_2));
  } while (abs(delta) > epsilon2 && --i > 0);
  return [x / (cos(phi) * (sqrt8 - 1 / sin(phi))), sign(y) * phi];
};

/* harmony default export */ var baker = (function() {
  return Object(d3_geo["geoProjection"])(bakerRaw)
      .scale(112.314);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/berghaus.js



function berghausRaw(lobes) {
  var k = 2 * pi / lobes;

  function forward(lambda, phi) {
    var p = Object(d3_geo["geoAzimuthalEquidistantRaw"])(lambda, phi);
    if (abs(lambda) > halfPi) { // back hemisphere
      var theta = atan2(p[1], p[0]),
          r = sqrt(p[0] * p[0] + p[1] * p[1]),
          theta0 = k * round((theta - halfPi) / k) + halfPi,
          α = atan2(sin(theta -= theta0), 2 - cos(theta)); // angle relative to lobe end
      theta = theta0 + asin(pi / r * sin(α)) - α;
      p[0] = r * cos(theta);
      p[1] = r * sin(theta);
    }
    return p;
  }

  forward.invert = function(x, y) {
    var r = sqrt(x * x + y * y);
    if (r > halfPi) {
      var theta = atan2(y, x),
          theta0 = k * round((theta - halfPi) / k) + halfPi,
          s = theta > theta0 ? -1 : 1,
          A = r * cos(theta0 - theta),
          cotα = 1 / tan(s * acos((A - pi) / sqrt(pi * (pi - 2 * A) + r * r)));
      theta = theta0 + 2 * atan((cotα + s * sqrt(cotα * cotα - 3)) / 3);
      x = r * cos(theta), y = r * sin(theta);
    }
    return d3_geo["geoAzimuthalEquidistantRaw"].invert(x, y);
  };

  return forward;
}

/* harmony default export */ var berghaus = (function() {
  var lobes = 5,
      m = Object(d3_geo["geoProjectionMutator"])(berghausRaw),
      p = m(lobes),
      projectionStream = p.stream,
      epsilon = 1e-2,
      cr = -cos(epsilon * radians),
      sr = sin(epsilon * radians);

  p.lobes = function(_) {
    return arguments.length ? m(lobes = +_) : lobes;
  };

  p.stream = function(stream) {
    var rotate = p.rotate(),
        rotateStream = projectionStream(stream),
        sphereStream = (p.rotate([0, 0]), projectionStream(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart(), sphereStream.lineStart();
      for (var i = 0, delta = 360 / lobes, delta0 = 2 * pi / lobes, phi = 90 - 180 / lobes, phi0 = halfPi; i < lobes; ++i, phi -= delta, phi0 -= delta0) {
        sphereStream.point(atan2(sr * cos(phi0), cr) * degrees, asin(sr * sin(phi0)) * degrees);
        if (phi < -90) {
          sphereStream.point(-90, -180 - phi - epsilon);
          sphereStream.point(-90, -180 - phi + epsilon);
        } else {
          sphereStream.point(90, phi + epsilon);
          sphereStream.point(90, phi - epsilon);
        }
      }
      sphereStream.lineEnd(), sphereStream.polygonEnd();
    };
    return rotateStream;
  };

  return p
      .scale(87.8076)
      .center([0, 17.1875])
      .clipAngle(180 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/mollweide.js



function mollweideBromleyTheta(cp, phi) {
  var cpsinPhi = cp * sin(phi), i = 30, delta;
  do phi -= delta = (phi + sin(phi) - cpsinPhi) / (1 + cos(phi));
  while (abs(delta) > math_epsilon && --i > 0);
  return phi / 2;
}

function mollweideBromleyRaw(cx, cy, cp) {

  function forward(lambda, phi) {
    return [cx * lambda * cos(phi = mollweideBromleyTheta(cp, phi)), cy * sin(phi)];
  }

  forward.invert = function(x, y) {
    return y = asin(y / cy), [x / (cx * cos(y)), asin((2 * y + sin(2 * y)) / cp)];
  };

  return forward;
}

var mollweideRaw = mollweideBromleyRaw(sqrt2 / halfPi, sqrt2, pi);

/* harmony default export */ var mollweide = (function() {
  return Object(d3_geo["geoProjection"])(mollweideRaw)
      .scale(169.529);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/boggs.js




var boggs_k = 2.00276,
    boggs_w = 1.11072;

function boggsRaw(lambda, phi) {
  var theta = mollweideBromleyTheta(pi, phi);
  return [boggs_k * lambda / (1 / cos(phi) + boggs_w / cos(theta)), (phi + sqrt2 * sin(theta)) / boggs_k];
}

boggsRaw.invert = function(x, y) {
  var ky = boggs_k * y, theta = y < 0 ? -quarterPi : quarterPi, i = 25, delta, phi;
  do {
    phi = ky - sqrt2 * sin(theta);
    theta -= delta = (sin(2 * theta) + 2 * theta - pi * sin(phi)) / (2 * cos(2 * theta) + 2 + pi * cos(phi) * sqrt2 * cos(theta));
  } while (abs(delta) > math_epsilon && --i > 0);
  phi = ky - sqrt2 * sin(theta);
  return [x * (1 / cos(phi) + boggs_w / cos(theta)) / boggs_k, phi];
};

/* harmony default export */ var boggs = (function() {
  return Object(d3_geo["geoProjection"])(boggsRaw)
      .scale(160.857);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/parallel1.js



/* harmony default export */ var parallel1 = (function(projectAt) {
  var phi0 = 0,
      m = Object(d3_geo["geoProjectionMutator"])(projectAt),
      p = m(phi0);

  p.parallel = function(_) {
    return arguments.length ? m(phi0 = _ * radians) : phi0 * degrees;
  };

  return p;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/sinusoidal.js



function sinusoidalRaw(lambda, phi) {
  return [lambda * cos(phi), phi];
}

sinusoidalRaw.invert = function(x, y) {
  return [x / cos(y), y];
};

/* harmony default export */ var sinusoidal = (function() {
  return Object(d3_geo["geoProjection"])(sinusoidalRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/bonne.js




function bonneRaw(phi0) {
  if (!phi0) return sinusoidalRaw;
  var cotPhi0 = 1 / tan(phi0);

  function forward(lambda, phi) {
    var rho = cotPhi0 + phi0 - phi,
        e = rho ? lambda * cos(phi) / rho : rho;
    return [rho * sin(e), cotPhi0 - rho * cos(e)];
  }

  forward.invert = function(x, y) {
    var rho = sqrt(x * x + (y = cotPhi0 - y) * y),
        phi = cotPhi0 + phi0 - rho;
    return [rho / cos(phi) * atan2(x, y), phi];
  };

  return forward;
}

/* harmony default export */ var bonne = (function() {
  return parallel1(bonneRaw)
      .scale(123.082)
      .center([0, 26.1441])
      .parallel(45);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/bottomley.js



function bottomleyRaw(sinPsi) {

  function forward(lambda, phi) {
    var rho = halfPi - phi,
        eta = rho ? lambda * sinPsi * sin(rho) / rho : rho;
    return [rho * sin(eta) / sinPsi, halfPi - rho * cos(eta)];
  }

  forward.invert = function(x, y) {
    var x1 = x * sinPsi,
        y1 = halfPi - y,
        rho = sqrt(x1 * x1 + y1 * y1),
        eta = atan2(x1, y1);
    return [(rho ? rho / sin(rho) : 1) * eta / sinPsi, halfPi - rho];
  };

  return forward;
}

/* harmony default export */ var bottomley = (function() {
  var sinPsi = 0.5,
      m = Object(d3_geo["geoProjectionMutator"])(bottomleyRaw),
      p = m(sinPsi);

  p.fraction = function(_) {
    return arguments.length ? m(sinPsi = +_) : sinPsi;
  };

  return p
      .scale(158.837);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/bromley.js




var bromleyRaw = mollweideBromleyRaw(1, 4 / pi, pi);

/* harmony default export */ var bromley = (function() {
  return Object(d3_geo["geoProjection"])(bromleyRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/chamberlin.js



// Azimuthal distance.
function chamberlin_distance(dPhi, c1, s1, c2, s2, dLambda) {
  var cosdLambda = cos(dLambda), r;
  if (abs(dPhi) > 1 || abs(dLambda) > 1) {
    r = acos(s1 * s2 + c1 * c2 * cosdLambda);
  } else {
    var sindPhi = sin(dPhi / 2), sindLambda = sin(dLambda / 2);
    r = 2 * asin(sqrt(sindPhi * sindPhi + c1 * c2 * sindLambda * sindLambda));
  }
  return abs(r) > math_epsilon ? [r, atan2(c2 * sin(dLambda), c1 * s2 - s1 * c2 * cosdLambda)] : [0, 0];
}

// Angle opposite a, and contained between sides of lengths b and c.
function angle(b, c, a) {
  return acos((b * b + c * c - a * a) / (2 * b * c));
}

// Normalize longitude.
function longitude(lambda) {
  return lambda - 2 * pi * floor((lambda + pi) / (2 * pi));
}

function chamberlinRaw(p0, p1, p2) {
  var points = [
    [p0[0], p0[1], sin(p0[1]), cos(p0[1])],
    [p1[0], p1[1], sin(p1[1]), cos(p1[1])],
    [p2[0], p2[1], sin(p2[1]), cos(p2[1])]
  ];

  for (var a = points[2], b, i = 0; i < 3; ++i, a = b) {
    b = points[i];
    a.v = chamberlin_distance(b[1] - a[1], a[3], a[2], b[3], b[2], b[0] - a[0]);
    a.point = [0, 0];
  }

  var beta0 = angle(points[0].v[0], points[2].v[0], points[1].v[0]),
      beta1 = angle(points[0].v[0], points[1].v[0], points[2].v[0]),
      beta2 = pi - beta0;

  points[2].point[1] = 0;
  points[0].point[0] = -(points[1].point[0] = points[0].v[0] / 2);

  var mean = [
    points[2].point[0] = points[0].point[0] + points[2].v[0] * cos(beta0),
    2 * (points[0].point[1] = points[1].point[1] = points[2].v[0] * sin(beta0))
  ];

  function forward(lambda, phi) {
    var sinPhi = sin(phi),
        cosPhi = cos(phi),
        v = new Array(3), i;

    // Compute distance and azimuth from control points.
    for (i = 0; i < 3; ++i) {
      var p = points[i];
      v[i] = chamberlin_distance(phi - p[1], p[3], p[2], cosPhi, sinPhi, lambda - p[0]);
      if (!v[i][0]) return p.point;
      v[i][1] = longitude(v[i][1] - p.v[1]);
    }

    // Arithmetic mean of interception points.
    var point = mean.slice();
    for (i = 0; i < 3; ++i) {
      var j = i == 2 ? 0 : i + 1;
      var a = angle(points[i].v[0], v[i][0], v[j][0]);
      if (v[i][1] < 0) a = -a;

      if (!i) {
        point[0] += v[i][0] * cos(a);
        point[1] -= v[i][0] * sin(a);
      } else if (i == 1) {
        a = beta1 - a;
        point[0] -= v[i][0] * cos(a);
        point[1] -= v[i][0] * sin(a);
      } else {
        a = beta2 - a;
        point[0] += v[i][0] * cos(a);
        point[1] += v[i][0] * sin(a);
      }
    }

    point[0] /= 3, point[1] /= 3;
    return point;
  }

  return forward;
}

function pointRadians(p) {
  return p[0] *= radians, p[1] *= radians, p;
}

function chamberlinAfrica() {
  return chamberlin([0, 22], [45, 22], [22.5, -22])
      .scale(380)
      .center([22.5, 2]);
}

function chamberlin(p0, p1, p2) { // TODO order matters!
  var c = Object(d3_geo["geoCentroid"])({type: "MultiPoint", coordinates: [p0, p1, p2]}),
      R = [-c[0], -c[1]],
      r = Object(d3_geo["geoRotation"])(R),
      p = Object(d3_geo["geoProjection"])(chamberlinRaw(pointRadians(r(p0)), pointRadians(r(p1)), pointRadians(r(p2)))).rotate(R),
      center = p.center;

  delete p.rotate;

  p.center = function(_) {
    return arguments.length ? center(r(_)) : r.invert(center());
  };

  return p
      .clipAngle(90);
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/collignon.js



function collignonRaw(lambda, phi) {
  var alpha = sqrt(1 - sin(phi));
  return [(2 / sqrtPi) * lambda * alpha, sqrtPi * (1 - alpha)];
}

collignonRaw.invert = function(x, y) {
  var lambda = (lambda = y / sqrtPi - 1) * lambda;
  return [lambda > 0 ? x * sqrt(pi / lambda) / 2 : 0, asin(1 - lambda)];
};

/* harmony default export */ var collignon = (function() {
  return Object(d3_geo["geoProjection"])(collignonRaw)
      .scale(95.6464)
      .center([0, 30]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/craig.js



function craigRaw(phi0) {
  var tanPhi0 = tan(phi0);

  function forward(lambda, phi) {
    return [lambda, (lambda ? lambda / sin(lambda) : 1) * (sin(phi) * cos(lambda) - tanPhi0 * cos(phi))];
  }

  forward.invert = tanPhi0 ? function(x, y) {
    if (x) y *= sin(x) / x;
    var cosλ = cos(x);
    return [x, 2 * atan2(sqrt(cosλ * cosλ + tanPhi0 * tanPhi0 - y * y) - cosλ, tanPhi0 - y)];
  } : function(x, y) {
    return [x, asin(x ? y * tan(x) / x : y)];
  };

  return forward;
}

/* harmony default export */ var craig = (function() {
  return parallel1(craigRaw)
      .scale(249.828)
      .clipAngle(90);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/craster.js



var sqrt3 = sqrt(3);

function crasterRaw(lambda, phi) {
  return [sqrt3 * lambda * (2 * cos(2 * phi / 3) - 1) / sqrtPi, sqrt3 * sqrtPi * sin(phi / 3)];
}

crasterRaw.invert = function(x, y) {
  var phi = 3 * asin(y / (sqrt3 * sqrtPi));
  return [sqrtPi * x / (sqrt3 * (2 * cos(2 * phi / 3) - 1)), phi];
};

/* harmony default export */ var craster = (function() {
  return Object(d3_geo["geoProjection"])(crasterRaw)
      .scale(156.19);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/cylindricalEqualArea.js



function cylindricalEqualAreaRaw(phi0) {
  var cosPhi0 = cos(phi0);

  function forward(lambda, phi) {
    return [lambda * cosPhi0, sin(phi) / cosPhi0];
  }

  forward.invert = function(x, y) {
    return [x / cosPhi0, asin(y * cosPhi0)];
  };

  return forward;
}

/* harmony default export */ var cylindricalEqualArea = (function() {
  return parallel1(cylindricalEqualAreaRaw)
      .parallel(38.58) // acos(sqrt(width / height / pi)) * radians
      .scale(195.044); // width / (sqrt(width / height / pi) * 2 * pi)
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/cylindricalStereographic.js



function cylindricalStereographicRaw(phi0) {
  var cosPhi0 = cos(phi0);

  function forward(lambda, phi) {
    return [lambda * cosPhi0, (1 + cosPhi0) * tan(phi / 2)];
  }

  forward.invert = function(x, y) {
    return [x / cosPhi0, atan(y / (1 + cosPhi0)) * 2];
  };

  return forward;
}

/* harmony default export */ var cylindricalStereographic = (function() {
  return parallel1(cylindricalStereographicRaw)
      .scale(124.75);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eckert1.js



function eckert1Raw(lambda, phi) {
  var alpha = sqrt(8 / (3 * pi));
  return [
    alpha * lambda * (1 - abs(phi) / pi),
    alpha * phi
  ];
}

eckert1Raw.invert = function(x, y) {
  var alpha = sqrt(8 / (3 * pi)),
      phi = y / alpha;
  return [
    x / (alpha * (1 - abs(phi) / pi)),
    phi
  ];
};

/* harmony default export */ var eckert1 = (function() {
  return Object(d3_geo["geoProjection"])(eckert1Raw)
      .scale(165.664);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eckert2.js



function eckert2Raw(lambda, phi) {
  var alpha = sqrt(4 - 3 * sin(abs(phi)));
  return [
    2 / sqrt(6 * pi) * lambda * alpha,
    sign(phi) * sqrt(2 * pi / 3) * (2 - alpha)
  ];
}

eckert2Raw.invert = function(x, y) {
  var alpha = 2 - abs(y) / sqrt(2 * pi / 3);
  return [
    x * sqrt(6 * pi) / (2 * alpha),
    sign(y) * asin((4 - alpha * alpha) / 3)
  ];
};

/* harmony default export */ var eckert2 = (function() {
  return Object(d3_geo["geoProjection"])(eckert2Raw)
      .scale(165.664);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eckert3.js



function eckert3Raw(lambda, phi) {
  var k = sqrt(pi * (4 + pi));
  return [
    2 / k * lambda * (1 + sqrt(1 - 4 * phi * phi / (pi * pi))),
    4 / k * phi
  ];
}

eckert3Raw.invert = function(x, y) {
  var k = sqrt(pi * (4 + pi)) / 2;
  return [
    x * k / (1 + sqrt(1 - y * y * (4 + pi) / (4 * pi))),
    y * k / 2
  ];
};

/* harmony default export */ var eckert3 = (function() {
  return Object(d3_geo["geoProjection"])(eckert3Raw)
      .scale(180.739);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eckert4.js



function eckert4Raw(lambda, phi) {
  var k = (2 + halfPi) * sin(phi);
  phi /= 2;
  for (var i = 0, delta = Infinity; i < 10 && abs(delta) > math_epsilon; i++) {
    var cosPhi = cos(phi);
    phi -= delta = (phi + sin(phi) * (cosPhi + 2) - k) / (2 * cosPhi * (1 + cosPhi));
  }
  return [
    2 / sqrt(pi * (4 + pi)) * lambda * (1 + cos(phi)),
    2 * sqrt(pi / (4 + pi)) * sin(phi)
  ];
}

eckert4Raw.invert = function(x, y) {
  var A = y * sqrt((4 + pi) / pi) / 2,
      k = asin(A),
      c = cos(k);
  return [
    x / (2 / sqrt(pi * (4 + pi)) * (1 + c)),
    asin((k + A * (c + 2)) / (2 + halfPi))
  ];
};

/* harmony default export */ var eckert4 = (function() {
  return Object(d3_geo["geoProjection"])(eckert4Raw)
      .scale(180.739);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eckert5.js



function eckert5Raw(lambda, phi) {
  return [
    lambda * (1 + cos(phi)) / sqrt(2 + pi),
    2 * phi / sqrt(2 + pi)
  ];
}

eckert5Raw.invert = function(x, y) {
  var k = sqrt(2 + pi),
      phi = y * k / 2;
  return [
    k * x / (1 + cos(phi)),
    phi
  ];
};

/* harmony default export */ var eckert5 = (function() {
  return Object(d3_geo["geoProjection"])(eckert5Raw)
      .scale(173.044);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eckert6.js



function eckert6Raw(lambda, phi) {
  var k = (1 + halfPi) * sin(phi);
  for (var i = 0, delta = Infinity; i < 10 && abs(delta) > math_epsilon; i++) {
    phi -= delta = (phi + sin(phi) - k) / (1 + cos(phi));
  }
  k = sqrt(2 + pi);
  return [
    lambda * (1 + cos(phi)) / k,
    2 * phi / k
  ];
}

eckert6Raw.invert = function(x, y) {
  var j = 1 + halfPi,
      k = sqrt(j / 2);
  return [
    x * 2 * k / (1 + cos(y *= k)),
    asin((y + sin(y)) / j)
  ];
};

/* harmony default export */ var eckert6 = (function() {
  return Object(d3_geo["geoProjection"])(eckert6Raw)
      .scale(173.044);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/eisenlohr.js




var eisenlohrK = 3 + 2 * sqrt2;

function eisenlohrRaw(lambda, phi) {
  var s0 = sin(lambda /= 2),
      c0 = cos(lambda),
      k = sqrt(cos(phi)),
      c1 = cos(phi /= 2),
      t = sin(phi) / (c1 + sqrt2 * c0 * k),
      c = sqrt(2 / (1 + t * t)),
      v = sqrt((sqrt2 * c1 + (c0 + s0) * k) / (sqrt2 * c1 + (c0 - s0) * k));
  return [
    eisenlohrK * (c * (v - 1 / v) - 2 * log(v)),
    eisenlohrK * (c * t * (v + 1 / v) - 2 * atan(t))
  ];
}

eisenlohrRaw.invert = function(x, y) {
  if (!(p = augustRaw.invert(x / 1.2, y * 1.065))) return null;
  var lambda = p[0], phi = p[1], i = 20, p;
  x /= eisenlohrK, y /= eisenlohrK;
  do {
    var _0 = lambda / 2,
        _1 = phi / 2,
        s0 = sin(_0),
        c0 = cos(_0),
        s1 = sin(_1),
        c1 = cos(_1),
        cos1 = cos(phi),
        k = sqrt(cos1),
        t = s1 / (c1 + sqrt2 * c0 * k),
        t2 = t * t,
        c = sqrt(2 / (1 + t2)),
        v0 = (sqrt2 * c1 + (c0 + s0) * k),
        v1 = (sqrt2 * c1 + (c0 - s0) * k),
        v2 = v0 / v1,
        v = sqrt(v2),
        vm1v = v - 1 / v,
        vp1v = v + 1 / v,
        fx = c * vm1v - 2 * log(v) - x,
        fy = c * t * vp1v - 2 * atan(t) - y,
        deltatDeltaLambda = s1 && sqrt1_2 * k * s0 * t2 / s1,
        deltatDeltaPhi = (sqrt2 * c0 * c1 + k) / (2 * (c1 + sqrt2 * c0 * k) * (c1 + sqrt2 * c0 * k) * k),
        deltacDeltat = -0.5 * t * c * c * c,
        deltacDeltaLambda = deltacDeltat * deltatDeltaLambda,
        deltacDeltaPhi = deltacDeltat * deltatDeltaPhi,
        A = (A = 2 * c1 + sqrt2 * k * (c0 - s0)) * A * v,
        deltavDeltaLambda = (sqrt2 * c0 * c1 * k + cos1) / A,
        deltavDeltaPhi = -(sqrt2 * s0 * s1) / (k * A),
        deltaxDeltaLambda = vm1v * deltacDeltaLambda - 2 * deltavDeltaLambda / v + c * (deltavDeltaLambda + deltavDeltaLambda / v2),
        deltaxDeltaPhi = vm1v * deltacDeltaPhi - 2 * deltavDeltaPhi / v + c * (deltavDeltaPhi + deltavDeltaPhi / v2),
        deltayDeltaLambda = t * vp1v * deltacDeltaLambda - 2 * deltatDeltaLambda / (1 + t2) + c * vp1v * deltatDeltaLambda + c * t * (deltavDeltaLambda - deltavDeltaLambda / v2),
        deltayDeltaPhi = t * vp1v * deltacDeltaPhi - 2 * deltatDeltaPhi / (1 + t2) + c * vp1v * deltatDeltaPhi + c * t * (deltavDeltaPhi - deltavDeltaPhi / v2),
        denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda;
    if (!denominator) break;
    var deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator,
        deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
    lambda -= deltaLambda;
    phi = max(-halfPi, min(halfPi, phi - deltaPhi));
  } while ((abs(deltaLambda) > math_epsilon || abs(deltaPhi) > math_epsilon) && --i > 0);
  return abs(abs(phi) - halfPi) < math_epsilon ? [0, phi] : i && [lambda, phi];
};

/* harmony default export */ var eisenlohr = (function() {
  return Object(d3_geo["geoProjection"])(eisenlohrRaw)
      .scale(62.5271);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/fahey.js



var faheyK = cos(35 * radians);

function faheyRaw(lambda, phi) {
  var t = tan(phi / 2);
  return [lambda * faheyK * sqrt(1 - t * t), (1 + faheyK) * t];
}

faheyRaw.invert = function(x, y) {
  var t = y / (1 + faheyK);
  return [x && x / (faheyK * sqrt(1 - t * t)), 2 * atan(t)];
};

/* harmony default export */ var fahey = (function() {
  return Object(d3_geo["geoProjection"])(faheyRaw)
      .scale(137.152);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/foucaut.js



function foucautRaw(lambda, phi) {
  var k = phi / 2, cosk = cos(k);
  return [ 2 * lambda / sqrtPi * cos(phi) * cosk * cosk, sqrtPi * tan(k)];
}

foucautRaw.invert = function(x, y) {
  var k = atan(y / sqrtPi), cosk = cos(k), phi = 2 * k;
  return [x * sqrtPi / 2 / (cos(phi) * cosk * cosk), phi];
};

/* harmony default export */ var foucaut = (function() {
  return Object(d3_geo["geoProjection"])(foucautRaw)
      .scale(135.264);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/gilbert.js



function gilbertForward(point) {
  return [point[0] / 2, asin(tan(point[1] / 2 * radians)) * degrees];
}

function gilbertInvert(point) {
  return [point[0] * 2, 2 * atan(sin(point[1] * radians)) * degrees];
}

/* harmony default export */ var src_gilbert = (function(projectionType) {
  if (projectionType == null) projectionType = d3_geo["geoOrthographic"];
  var projection = projectionType(),
      equirectangular = Object(d3_geo["geoEquirectangular"])().scale(degrees).precision(0).clipAngle(null).translate([0, 0]); // antimeridian cutting

  function gilbert(point) {
    return projection(gilbertForward(point));
  }

  if (projection.invert) gilbert.invert = function(point) {
    return gilbertInvert(projection.invert(point));
  };

  gilbert.stream = function(stream) {
    var s1 = projection.stream(stream), s0 = equirectangular.stream({
      point: function(lambda, phi) { s1.point(lambda / 2, asin(tan(-phi / 2 * radians)) * degrees); },
      lineStart: function() { s1.lineStart(); },
      lineEnd: function() { s1.lineEnd(); },
      polygonStart: function() { s1.polygonStart(); },
      polygonEnd: function() { s1.polygonEnd(); }
    });
    s0.sphere = s1.sphere;
    return s0;
  };

  function property(name) {
    gilbert[name] = function(_) {
      return arguments.length ? (projection[name](_), gilbert) : projection[name]();
    };
  }

  gilbert.rotate = function(_) {
    return arguments.length ? (equirectangular.rotate(_), gilbert) : equirectangular.rotate();
  };

  gilbert.center = function(_) {
    return arguments.length ? (projection.center(gilbertForward(_)), gilbert) : gilbertInvert(projection.center());
  };

  property("clipAngle");
  property("clipExtent");
  property("scale");
  property("translate");
  property("precision");

  return gilbert
      .scale(249.5);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/gingery.js



function gingeryRaw(rho, n) {
  var k = 2 * pi / n,
      rho2 = rho * rho;

  function forward(lambda, phi) {
    var p = Object(d3_geo["geoAzimuthalEquidistantRaw"])(lambda, phi),
        x = p[0],
        y = p[1],
        r2 = x * x + y * y;

    if (r2 > rho2) {
      var r = sqrt(r2),
          theta = atan2(y, x),
          theta0 = k * round(theta / k),
          alpha = theta - theta0,
          rhoCosAlpha = rho * cos(alpha),
          k_ = (rho * sin(alpha) - alpha * sin(rhoCosAlpha)) / (halfPi - rhoCosAlpha),
          s_ = gingeryLength(alpha, k_),
          e = (pi - rho) / gingeryIntegrate(s_, rhoCosAlpha, pi);

      x = r;
      var i = 50, delta;
      do {
        x -= delta = (rho + gingeryIntegrate(s_, rhoCosAlpha, x) * e - r) / (s_(x) * e);
      } while (abs(delta) > math_epsilon && --i > 0);

      y = alpha * sin(x);
      if (x < halfPi) y -= k_ * (x - halfPi);

      var s = sin(theta0),
          c = cos(theta0);
      p[0] = x * c - y * s;
      p[1] = x * s + y * c;
    }
    return p;
  }

  forward.invert = function(x, y) {
    var r2 = x * x + y * y;
    if (r2 > rho2) {
      var r = sqrt(r2),
          theta = atan2(y, x),
          theta0 = k * round(theta / k),
          dTheta = theta - theta0;

      x = r * cos(dTheta);
      y = r * sin(dTheta);

      var x_halfPi = x - halfPi,
          sinx = sin(x),
          alpha = y / sinx,
          delta = x < halfPi ? Infinity : 0,
          i = 10;

      while (true) {
        var rhosinAlpha = rho * sin(alpha),
            rhoCosAlpha = rho * cos(alpha),
            sinRhoCosAlpha = sin(rhoCosAlpha),
            halfPi_RhoCosAlpha = halfPi - rhoCosAlpha,
            k_ = (rhosinAlpha - alpha * sinRhoCosAlpha) / halfPi_RhoCosAlpha,
            s_ = gingeryLength(alpha, k_);

        if (abs(delta) < epsilon2 || !--i) break;

        alpha -= delta = (alpha * sinx - k_ * x_halfPi - y) / (
          sinx - x_halfPi * 2 * (
            halfPi_RhoCosAlpha * (rhoCosAlpha + alpha * rhosinAlpha * cos(rhoCosAlpha) - sinRhoCosAlpha) -
            rhosinAlpha * (rhosinAlpha - alpha * sinRhoCosAlpha)
          ) / (halfPi_RhoCosAlpha * halfPi_RhoCosAlpha));
      }
      r = rho + gingeryIntegrate(s_, rhoCosAlpha, x) * (pi - rho) / gingeryIntegrate(s_, rhoCosAlpha, pi);
      theta = theta0 + alpha;
      x = r * cos(theta);
      y = r * sin(theta);
    }
    return d3_geo["geoAzimuthalEquidistantRaw"].invert(x, y);
  };

  return forward;
}

function gingeryLength(alpha, k) {
  return function(x) {
    var y_ = alpha * cos(x);
    if (x < halfPi) y_ -= k;
    return sqrt(1 + y_ * y_);
  };
}

// Numerical integration: trapezoidal rule.
function gingeryIntegrate(f, a, b) {
  var n = 50,
      h = (b - a) / n,
      s = f(a) + f(b);
  for (var i = 1, x = a; i < n; ++i) s += 2 * f(x += h);
  return s * 0.5 * h;
}

/* harmony default export */ var gingery = (function() {
  var n = 6,
      rho = 30 * radians,
      cRho = cos(rho),
      sRho = sin(rho),
      m = Object(d3_geo["geoProjectionMutator"])(gingeryRaw),
      p = m(rho, n),
      stream_ = p.stream,
      epsilon = 1e-2,
      cr = -cos(epsilon * radians),
      sr = sin(epsilon * radians);

  p.radius = function(_) {
    if (!arguments.length) return rho * degrees;
    cRho = cos(rho = _ * radians);
    sRho = sin(rho);
    return m(rho, n);
  };

  p.lobes = function(_) {
    if (!arguments.length) return n;
    return m(rho, n = +_);
  };

  p.stream = function(stream) {
    var rotate = p.rotate(),
        rotateStream = stream_(stream),
        sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart(), sphereStream.lineStart();
      for (var i = 0, delta = 2 * pi / n, phi = 0; i < n; ++i, phi -= delta) {
        sphereStream.point(atan2(sr * cos(phi), cr) * degrees, asin(sr * sin(phi)) * degrees);
        sphereStream.point(atan2(sRho * cos(phi - delta / 2), cRho) * degrees, asin(sRho * sin(phi - delta / 2)) * degrees);
      }
      sphereStream.lineEnd(), sphereStream.polygonEnd();
    };
    return rotateStream;
  };

  return p
      .rotate([90, -40])
      .scale(91.7095)
      .clipAngle(180 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/ginzburgPolyconic.js


/* harmony default export */ var ginzburgPolyconic = (function(a, b, c, d, e, f, g, h) {
  if (arguments.length < 8) h = 0;

  function forward(lambda, phi) {
    if (!phi) return [a * lambda / pi, 0];
    var phi2 = phi * phi,
        xB = a + phi2 * (b + phi2 * (c + phi2 * d)),
        yB = phi * (e - 1 + phi2 * (f - h + phi2 * g)),
        m = (xB * xB + yB * yB) / (2 * yB),
        alpha = lambda * asin(xB / m) / pi;
    return [m * sin(alpha), phi * (1 + phi2 * h) + m * (1 - cos(alpha))];
  }

  forward.invert = function(x, y) {
    var lambda = pi * x / a,
        phi = y,
        deltaLambda, deltaPhi, i = 50;
    do {
      var phi2 = phi * phi,
          xB = a + phi2 * (b + phi2 * (c + phi2 * d)),
          yB = phi * (e - 1 + phi2 * (f - h + phi2 * g)),
          p = xB * xB + yB * yB,
          q = 2 * yB,
          m = p / q,
          m2 = m * m,
          dAlphadLambda = asin(xB / m) / pi,
          alpha = lambda * dAlphadLambda,
          xB2 = xB * xB,
          dxBdPhi = (2 * b + phi2 * (4 * c + phi2 * 6 * d)) * phi,
          dyBdPhi = e + phi2 * (3 * f + phi2 * 5 * g),
          dpdPhi = 2 * (xB * dxBdPhi + yB * (dyBdPhi - 1)),
          dqdPhi = 2 * (dyBdPhi - 1),
          dmdPhi = (dpdPhi * q - p * dqdPhi) / (q * q),
          cosAlpha = cos(alpha),
          sinAlpha = sin(alpha),
          mcosAlpha = m * cosAlpha,
          msinAlpha = m * sinAlpha,
          dAlphadPhi = ((lambda / pi) * (1 / sqrt(1 - xB2 / m2)) * (dxBdPhi * m - xB * dmdPhi)) / m2,
          fx = msinAlpha - x,
          fy = phi * (1 + phi2 * h) + m - mcosAlpha - y,
          deltaxDeltaPhi = dmdPhi * sinAlpha + mcosAlpha * dAlphadPhi,
          deltaxDeltaLambda = mcosAlpha * dAlphadLambda,
          deltayDeltaPhi = 1 + dmdPhi - (dmdPhi * cosAlpha - msinAlpha * dAlphadPhi),
          deltayDeltaLambda = msinAlpha * dAlphadLambda,
          denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda;
      if (!denominator) break;
      lambda -= deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator;
      phi -= deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
    } while ((abs(deltaLambda) > math_epsilon || abs(deltaPhi) > math_epsilon) && --i > 0);
    return [lambda, phi];
  };

  return forward;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/ginzburg4.js



var ginzburg4Raw = ginzburgPolyconic(2.8284, -1.6988, 0.75432, -0.18071, 1.76003, -0.38914, 0.042555);

/* harmony default export */ var ginzburg4 = (function() {
  return Object(d3_geo["geoProjection"])(ginzburg4Raw)
      .scale(149.995);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/ginzburg5.js



var ginzburg5Raw = ginzburgPolyconic(2.583819, -0.835827, 0.170354, -0.038094, 1.543313, -0.411435,0.082742);

/* harmony default export */ var ginzburg5 = (function() {
  return Object(d3_geo["geoProjection"])(ginzburg5Raw)
      .scale(153.93);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/ginzburg6.js




var ginzburg6Raw = ginzburgPolyconic(5 / 6 * pi, -0.62636, -0.0344, 0, 1.3493, -0.05524, 0, 0.045);

/* harmony default export */ var ginzburg6 = (function() {
  return Object(d3_geo["geoProjection"])(ginzburg6Raw)
      .scale(130.945);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/ginzburg8.js



function ginzburg8Raw(lambda, phi) {
  var lambda2 = lambda * lambda,
      phi2 = phi * phi;
  return [
    lambda * (1 - 0.162388 * phi2) * (0.87 - 0.000952426 * lambda2 * lambda2),
    phi * (1 + phi2 / 12)
  ];
}

ginzburg8Raw.invert = function(x, y) {
  var lambda = x,
      phi = y,
      i = 50, delta;
  do {
    var phi2 = phi * phi;
    phi -= delta = (phi * (1 + phi2 / 12) - y) / (1 + phi2 / 4);
  } while (abs(delta) > math_epsilon && --i > 0);
  i = 50;
  x /= 1 -0.162388 * phi2;
  do {
    var lambda4 = (lambda4 = lambda * lambda) * lambda4;
    lambda -= delta = (lambda * (0.87 - 0.000952426 * lambda4) - x) / (0.87 - 0.00476213 * lambda4);
  } while (abs(delta) > math_epsilon && --i > 0);
  return [lambda, phi];
};

/* harmony default export */ var ginzburg8 = (function() {
  return Object(d3_geo["geoProjection"])(ginzburg8Raw)
      .scale(131.747);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/ginzburg9.js



var ginzburg9Raw = ginzburgPolyconic(2.6516, -0.76534, 0.19123, -0.047094, 1.36289, -0.13965,0.031762);

/* harmony default export */ var ginzburg9 = (function() {
  return Object(d3_geo["geoProjection"])(ginzburg9Raw)
      .scale(131.087);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/square.js


/* harmony default export */ var square = (function(project) {
  var dx = project(halfPi, 0)[0] - project(-halfPi, 0)[0];

  function projectSquare(lambda, phi) {
    var s = lambda > 0 ? -0.5 : 0.5,
        point = project(lambda + s * pi, phi);
    point[0] -= s * dx;
    return point;
  }

  if (project.invert) projectSquare.invert = function(x, y) {
    var s = x > 0 ? -0.5 : 0.5,
        location = project.invert(x + s * dx, y),
        lambda = location[0] - s * pi;
    if (lambda < -pi) lambda += 2 * pi;
    else if (lambda > pi) lambda -= 2 * pi;
    location[0] = lambda;
    return location;
  };

  return projectSquare;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/gringorten.js




function gringortenRaw(lambda, phi) {
  var sLambda = sign(lambda),
      sPhi = sign(phi),
      cosPhi = cos(phi),
      x = cos(lambda) * cosPhi,
      y = sin(lambda) * cosPhi,
      z = sin(sPhi * phi);
  lambda = abs(atan2(y, z));
  phi = asin(x);
  if (abs(lambda - halfPi) > math_epsilon) lambda %= halfPi;
  var point = gringortenHexadecant(lambda > pi / 4 ? halfPi - lambda : lambda, phi);
  if (lambda > pi / 4) z = point[0], point[0] = -point[1], point[1] = -z;
  return (point[0] *= sLambda, point[1] *= -sPhi, point);
}

gringortenRaw.invert = function(x, y) {
  var sx = sign(x),
      sy = sign(y),
      x0 = -sx * x,
      y0 = -sy * y,
      t = y0 / x0 < 1,
      p = gringortenHexadecantInvert(t ? y0 : x0, t ? x0 : y0),
      lambda = p[0],
      phi = p[1],
      cosPhi = cos(phi);
  if (t) lambda = -halfPi - lambda;
  return [sx * (atan2(sin(lambda) * cosPhi, -sin(phi)) + pi), sy * asin(cos(lambda) * cosPhi)];
};

function gringortenHexadecant(lambda, phi) {
  if (phi === halfPi) return [0, 0];

  var sinPhi = sin(phi),
      r = sinPhi * sinPhi,
      r2 = r * r,
      j = 1 + r2,
      k = 1 + 3 * r2,
      q = 1 - r2,
      z = asin(1 / sqrt(j)),
      v = q + r * j * z,
      p2 = (1 - sinPhi) / v,
      p = sqrt(p2),
      a2 = p2 * j,
      a = sqrt(a2),
      h = p * q,
      x,
      i;

  if (lambda === 0) return [0, -(h + r * a)];

  var cosPhi = cos(phi),
      secPhi = 1 / cosPhi,
      drdPhi = 2 * sinPhi * cosPhi,
      dvdPhi = (-3 * r + z * k) * drdPhi,
      dp2dPhi = (-v * cosPhi - (1 - sinPhi) * dvdPhi) / (v * v),
      dpdPhi = (0.5 * dp2dPhi) / p,
      dhdPhi = q * dpdPhi - 2 * r * p * drdPhi,
      dra2dPhi = r * j * dp2dPhi + p2 * k * drdPhi,
      mu = -secPhi * drdPhi,
      nu = -secPhi * dra2dPhi,
      zeta = -2 * secPhi * dhdPhi,
      lambda1 = 4 * lambda / pi,
      delta;

  // Slower but accurate bisection method.
  if (lambda > 0.222 * pi || phi < pi / 4 && lambda > 0.175 * pi) {
    x = (h + r * sqrt(a2 * (1 + r2) - h * h)) / (1 + r2);
    if (lambda > pi / 4) return [x, x];
    var x1 = x, x0 = 0.5 * x;
    x = 0.5 * (x0 + x1), i = 50;
    do {
      var g = sqrt(a2 - x * x),
          f = (x * (zeta + mu * g) + nu * asin(x / a)) - lambda1;
      if (!f) break;
      if (f < 0) x0 = x;
      else x1 = x;
      x = 0.5 * (x0 + x1);
    } while (abs(x1 - x0) > math_epsilon && --i > 0);
  }

  // Newton-Raphson.
  else {
    x = math_epsilon, i = 25;
    do {
      var x2 = x * x,
          g2 = sqrt(a2 - x2),
          zetaMug = zeta + mu * g2,
          f2 = x * zetaMug + nu * asin(x / a) - lambda1,
          df = zetaMug + (nu - mu * x2) / g2;
      x -= delta = g2 ? f2 / df : 0;
    } while (abs(delta) > math_epsilon && --i > 0);
  }

  return [x, -h - r * sqrt(a2 - x * x)];
}

function gringortenHexadecantInvert(x, y) {
  var x0 = 0,
      x1 = 1,
      r = 0.5,
      i = 50;

  while (true) {
    var r2 = r * r,
        sinPhi = sqrt(r),
        z = asin(1 / sqrt(1 + r2)),
        v = (1 - r2) + r * (1 + r2) * z,
        p2 = (1 - sinPhi) / v,
        p = sqrt(p2),
        a2 = p2 * (1 + r2),
        h = p * (1 - r2),
        g2 = a2 - x * x,
        g = sqrt(g2),
        y0 = y + h + r * g;
    if (abs(x1 - x0) < epsilon2 || --i === 0 || y0 === 0) break;
    if (y0 > 0) x0 = r;
    else x1 = r;
    r = 0.5 * (x0 + x1);
  }

  if (!i) return null;

  var phi = asin(sinPhi),
      cosPhi = cos(phi),
      secPhi = 1 / cosPhi,
      drdPhi = 2 * sinPhi * cosPhi,
      dvdPhi = (-3 * r + z * (1 + 3 * r2)) * drdPhi,
      dp2dPhi = (-v * cosPhi - (1 - sinPhi) * dvdPhi) / (v * v),
      dpdPhi = 0.5 * dp2dPhi / p,
      dhdPhi = (1 - r2) * dpdPhi - 2 * r * p * drdPhi,
      zeta = -2 * secPhi * dhdPhi,
      mu = -secPhi * drdPhi,
      nu = -secPhi * (r * (1 + r2) * dp2dPhi + p2 * (1 + 3 * r2) * drdPhi);

  return [pi / 4 * (x * (zeta + mu * g) + nu * asin(x / sqrt(a2))), phi];
}

/* harmony default export */ var gringorten = (function() {
  return Object(d3_geo["geoProjection"])(square(gringortenRaw))
      .scale(239.75);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/elliptic.js


// Returns [sn, cn, dn](u + iv|m).
function ellipticJi(u, v, m) {
  var a, b, c;
  if (!u) {
    b = ellipticJ(v, 1 - m);
    return [
      [0, b[0] / b[1]],
      [1 / b[1], 0],
      [b[2] / b[1], 0]
    ];
  }
  a = ellipticJ(u, m);
  if (!v) return [[a[0], 0], [a[1], 0], [a[2], 0]];
  b = ellipticJ(v, 1 - m);
  c = b[1] * b[1] + m * a[0] * a[0] * b[0] * b[0];
  return [
    [a[0] * b[2] / c, a[1] * a[2] * b[0] * b[1] / c],
    [a[1] * b[1] / c, -a[0] * a[2] * b[0] * b[2] / c],
    [a[2] * b[1] * b[2] / c, -m * a[0] * a[1] * b[0] / c]
  ];
}

// Returns [sn, cn, dn, ph](u|m).
function ellipticJ(u, m) {
  var ai, b, phi, t, twon;
  if (m < math_epsilon) {
    t = sin(u);
    b = cos(u);
    ai = m * (u - t * b) / 4;
    return [
      t - ai * b,
      b + ai * t,
      1 - m * t * t / 2,
      u - ai
    ];
  }
  if (m >= 1 - math_epsilon) {
    ai = (1 - m) / 4;
    b = cosh(u);
    t = tanh(u);
    phi = 1 / b;
    twon = b * sinh(u);
    return [
      t + ai * (twon - u) / (b * b),
      phi - ai * t * phi * (twon - u),
      phi + ai * t * phi * (twon + u),
      2 * atan(exp(u)) - halfPi + ai * (twon - u) / b
    ];
  }

  var a = [1, 0, 0, 0, 0, 0, 0, 0, 0],
      c = [sqrt(m), 0, 0, 0, 0, 0, 0, 0, 0],
      i = 0;
  b = sqrt(1 - m);
  twon = 1;

  while (abs(c[i] / a[i]) > math_epsilon && i < 8) {
    ai = a[i++];
    c[i] = (ai - b) / 2;
    a[i] = (ai + b) / 2;
    b = sqrt(ai * b);
    twon *= 2;
  }

  phi = twon * a[i] * u;
  do {
    t = c[i] * sin(b = phi) / a[i];
    phi = (asin(t) + phi) / 2;
  } while (--i);

  return [sin(phi), t = cos(phi), t / cos(phi - b), phi];
}

// Calculate F(phi+iPsi|m).
// See Abramowitz and Stegun, 17.4.11.
function ellipticFi(phi, psi, m) {
  var r = abs(phi),
      i = abs(psi),
      sinhPsi = sinh(i);
  if (r) {
    var cscPhi = 1 / sin(r),
        cotPhi2 = 1 / (tan(r) * tan(r)),
        b = -(cotPhi2 + m * (sinhPsi * sinhPsi * cscPhi * cscPhi) - 1 + m),
        c = (m - 1) * cotPhi2,
        cotLambda2 = (-b + sqrt(b * b - 4 * c)) / 2;
    return [
      ellipticF(atan(1 / sqrt(cotLambda2)), m) * sign(phi),
      ellipticF(atan(sqrt((cotLambda2 / cotPhi2 - 1) / m)), 1 - m) * sign(psi)
    ];
  }
  return [
    0,
    ellipticF(atan(sinhPsi), 1 - m) * sign(psi)
  ];
}

// Calculate F(phi|m) where m = k² = sin²α.
// See Abramowitz and Stegun, 17.6.7.
function ellipticF(phi, m) {
  if (!m) return phi;
  if (m === 1) return log(tan(phi / 2 + quarterPi));
  var a = 1,
      b = sqrt(1 - m),
      c = sqrt(m);
  for (var i = 0; abs(c) > math_epsilon; i++) {
    if (phi % pi) {
      var dPhi = atan(b * tan(phi) / a);
      if (dPhi < 0) dPhi += pi;
      phi += dPhi + ~~(phi / pi) * pi;
    } else phi += phi;
    c = (a + b) / 2;
    b = sqrt(a * b);
    c = ((a = c) - b) / 2;
  }
  return phi / (pow(2, i) * a);
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/guyou.js





function guyouRaw(lambda, phi) {
  var k_ = (sqrt2 - 1) / (sqrt2 + 1),
      k = sqrt(1 - k_ * k_),
      K = ellipticF(halfPi, k * k),
      f = -1,
      psi = log(tan(pi / 4 + abs(phi) / 2)),
      r = exp(f * psi) / sqrt(k_),
      at = guyouComplexAtan(r * cos(f * lambda), r * sin(f * lambda)),
      t = ellipticFi(at[0], at[1], k * k);
  return [-t[1], (phi >= 0 ? 1 : -1) * (0.5 * K - t[0])];
}

function guyouComplexAtan(x, y) {
  var x2 = x * x,
      y_1 = y + 1,
      t = 1 - x2 - y * y;
  return [
   0.5 * ((x >= 0 ? halfPi : -halfPi) - atan2(t, 2 * x)),
    -0.25 * log(t * t + 4 * x2) +0.5 * log(y_1 * y_1 + x2)
  ];
}

function guyouComplexDivide(a, b) {
  var denominator = b[0] * b[0] + b[1] * b[1];
  return [
    (a[0] * b[0] + a[1] * b[1]) / denominator,
    (a[1] * b[0] - a[0] * b[1]) / denominator
  ];
}

guyouRaw.invert = function(x, y) {
  var k_ = (sqrt2 - 1) / (sqrt2 + 1),
      k = sqrt(1 - k_ * k_),
      K = ellipticF(halfPi, k * k),
      f = -1,
      j = ellipticJi(0.5 * K - y, -x, k * k),
      tn = guyouComplexDivide(j[0], j[1]),
      lambda = atan2(tn[1], tn[0]) / f;
  return [
    lambda,
    2 * atan(exp(0.5 / f * log(k_ * tn[0] * tn[0] + k_ * tn[1] * tn[1]))) - halfPi
  ];
};

/* harmony default export */ var guyou = (function() {
  return Object(d3_geo["geoProjection"])(square(guyouRaw))
      .scale(151.496);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/hammer.js



function hammerRaw(A, B) {
  if (arguments.length < 2) B = A;
  if (B === 1) return d3_geo["geoAzimuthalEqualAreaRaw"];
  if (B === Infinity) return hammerQuarticAuthalicRaw;

  function forward(lambda, phi) {
    var coordinates = Object(d3_geo["geoAzimuthalEqualAreaRaw"])(lambda / B, phi);
    coordinates[0] *= A;
    return coordinates;
  }

  forward.invert = function(x, y) {
    var coordinates = d3_geo["geoAzimuthalEqualAreaRaw"].invert(x / A, y);
    coordinates[0] *= B;
    return coordinates;
  };

  return forward;
}

function hammerQuarticAuthalicRaw(lambda, phi) {
  return [
    lambda * cos(phi) / cos(phi /= 2),
    2 * sin(phi)
  ];
}

hammerQuarticAuthalicRaw.invert = function(x, y) {
  var phi = 2 * asin(y / 2);
  return [
    x * cos(phi / 2) / cos(phi),
    phi
  ];
};

/* harmony default export */ var hammer = (function() {
  var B = 2,
      m = Object(d3_geo["geoProjectionMutator"])(hammerRaw),
      p = m(B);

  p.coefficient = function(_) {
    if (!arguments.length) return B;
    return m(B = +_);
  };

  return p
    .scale(169.529);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/hammerRetroazimuthal.js



function hammerRetroazimuthalRaw(phi0) {
  var sinPhi0 = sin(phi0),
      cosPhi0 = cos(phi0),
      rotate = hammerRetroazimuthalRotation(phi0);

  rotate.invert = hammerRetroazimuthalRotation(-phi0);

  function forward(lambda, phi) {
    var p = rotate(lambda, phi);
    lambda = p[0], phi = p[1];
    var sinPhi = sin(phi),
        cosPhi = cos(phi),
        cosLambda = cos(lambda),
        z = acos(sinPhi0 * sinPhi + cosPhi0 * cosPhi * cosLambda),
        sinz = sin(z),
        K = abs(sinz) > math_epsilon ? z / sinz : 1;
    return [
      K * cosPhi0 * sin(lambda),
      (abs(lambda) > halfPi ? K : -K) // rotate for back hemisphere
        * (sinPhi0 * cosPhi - cosPhi0 * sinPhi * cosLambda)
    ];
  }

  forward.invert = function(x, y) {
    var rho = sqrt(x * x + y * y),
        sinz = -sin(rho),
        cosz = cos(rho),
        a = rho * cosz,
        b = -y * sinz,
        c = rho * sinPhi0,
        d = sqrt(a * a + b * b - c * c),
        phi = atan2(a * c + b * d, b * c - a * d),
        lambda = (rho > halfPi ? -1 : 1) * atan2(x * sinz, rho * cos(phi) * cosz + y * sin(phi) * sinz);
    return rotate.invert(lambda, phi);
  };

  return forward;
}

// Latitudinal rotation by phi0.
// Temporary hack until D3 supports arbitrary small-circle clipping origins.
function hammerRetroazimuthalRotation(phi0) {
  var sinPhi0 = sin(phi0),
      cosPhi0 = cos(phi0);

  return function(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi);
    return [
      atan2(y, x * cosPhi0 - z * sinPhi0),
      asin(z * cosPhi0 + x * sinPhi0)
    ];
  };
}

/* harmony default export */ var hammerRetroazimuthal = (function() {
  var phi0 = 0,
      m = Object(d3_geo["geoProjectionMutator"])(hammerRetroazimuthalRaw),
      p = m(phi0),
      rotate_ = p.rotate,
      stream_ = p.stream,
      circle = Object(d3_geo["geoCircle"])();

  p.parallel = function(_) {
    if (!arguments.length) return phi0 * degrees;
    var r = p.rotate();
    return m(phi0 = _ * radians).rotate(r);
  };

  // Temporary hack; see hammerRetroazimuthalRotation.
  p.rotate = function(_) {
    if (!arguments.length) return (_ = rotate_.call(p), _[1] += phi0 * degrees, _);
    rotate_.call(p, [_[0], _[1] - phi0 * degrees]);
    circle.center([-_[0], -_[1]]);
    return p;
  };

  p.stream = function(stream) {
    stream = stream_(stream);
    stream.sphere = function() {
      stream.polygonStart();
      var epsilon = 1e-2,
          ring = circle.radius(90 - epsilon)().coordinates[0],
          n = ring.length - 1,
          i = -1,
          p;
      stream.lineStart();
      while (++i < n) stream.point((p = ring[i])[0], p[1]);
      stream.lineEnd();
      ring = circle.radius(90 + epsilon)().coordinates[0];
      n = ring.length - 1;
      stream.lineStart();
      while (--i >= 0) stream.point((p = ring[i])[0], p[1]);
      stream.lineEnd();
      stream.polygonEnd();
    };
    return stream;
  };

  return p
      .scale(79.4187)
      .parallel(45)
      .clipAngle(180 - 1e-3);
});

// EXTERNAL MODULE: ./node_modules/d3-array/src/index.js + 31 modules
var src = __webpack_require__(2);

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/healpix.js






var healpixParallel = 41 + 48 / 36 + 37 / 3600, // for K=3; TODO automate
    healpixLambert = cylindricalEqualAreaRaw(0);

function healpixRaw(H) {
  var phi0 = healpixParallel * radians,
      dx = collignonRaw(pi, phi0)[0] - collignonRaw(-pi, phi0)[0],
      y0 = healpixLambert(0, phi0)[1],
      y1 = collignonRaw(0, phi0)[1],
      dy1 = sqrtPi - y1,
      k = tau / H,
      w = 4 / tau,
      h = y0 + (dy1 * dy1 * 4) / tau;

  function forward(lambda, phi) {
    var point,
        phi2 = abs(phi);
    if (phi2 > phi0) {
      var i = min(H - 1, max(0, floor((lambda + pi) / k)));
      lambda += pi * (H - 1) / H - i * k;
      point = collignonRaw(lambda, phi2);
      point[0] = point[0] * tau / dx - tau * (H - 1) / (2 * H) + i * tau / H;
      point[1] = y0 + (point[1] - y1) * 4 * dy1 / tau;
      if (phi < 0) point[1] = -point[1];
    } else {
      point = healpixLambert(lambda, phi);
    }
    point[0] *= w, point[1] /= h;
    return point;
  }

  forward.invert = function(x, y) {
    x /= w, y *= h;
    var y2 = abs(y);
    if (y2 > y0) {
      var i = min(H - 1, max(0, floor((x + pi) / k)));
      x = (x + pi * (H - 1) / H - i * k) * dx / tau;
      var point = collignonRaw.invert(x, 0.25 * (y2 - y0) * tau / dy1 + y1);
      point[0] -= pi * (H - 1) / H - i * k;
      if (y < 0) point[1] = -point[1];
      return point;
    }
    return healpixLambert.invert(x, y);
  };

  return forward;
}

function healpix_sphere(step) {
  return {
    type: "Polygon",
    coordinates: [
      Object(src["b" /* range */])(-180, 180 + step / 2, step).map(function(x, i) { return [x, i & 1 ? 90 - 1e-6 : healpixParallel]; })
      .concat(Object(src["b" /* range */])(180, -180 - step / 2, -step).map(function(x, i) { return [x, i & 1 ? -90 + 1e-6 : -healpixParallel]; }))
    ]
  };
}

/* harmony default export */ var healpix = (function() {
  var H = 4,
      m = Object(d3_geo["geoProjectionMutator"])(healpixRaw),
      p = m(H),
      stream_ = p.stream;

  p.lobes = function(_) {
    return arguments.length ? m(H = +_) : H;
  };

  p.stream = function(stream) {
    var rotate = p.rotate(),
        rotateStream = stream_(stream),
        sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() { Object(d3_geo["geoStream"])(healpix_sphere(180 / H), sphereStream); };
    return rotateStream;
  };

  return p
      .scale(239.75);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/hill.js



function hillRaw(K) {
  var L = 1 + K,
      sinBt = sin(1 / L),
      Bt = asin(sinBt),
      A = 2 * sqrt(pi / (B = pi + 4 * Bt * L)),
      B,
      rho0 = 0.5 * A * (L + sqrt(K * (2 + K))),
      K2 = K * K,
      L2 = L * L;

  function forward(lambda, phi) {
    var t = 1 - sin(phi),
        rho,
        omega;
    if (t && t < 2) {
      var theta = halfPi - phi, i = 25, delta;
      do {
        var sinTheta = sin(theta),
            cosTheta = cos(theta),
            Bt_Bt1 = Bt + atan2(sinTheta, L - cosTheta),
            C = 1 + L2 - 2 * L * cosTheta;
        theta -= delta = (theta - K2 * Bt - L * sinTheta + C * Bt_Bt1 -0.5 * t * B) / (2 * L * sinTheta * Bt_Bt1);
      } while (abs(delta) > epsilon2 && --i > 0);
      rho = A * sqrt(C);
      omega = lambda * Bt_Bt1 / pi;
    } else {
      rho = A * (K + t);
      omega = lambda * Bt / pi;
    }
    return [
      rho * sin(omega),
      rho0 - rho * cos(omega)
    ];
  }

  forward.invert = function(x, y) {
    var rho2 = x * x + (y -= rho0) * y,
        cosTheta = (1 + L2 - rho2 / (A * A)) / (2 * L),
        theta = acos(cosTheta),
        sinTheta = sin(theta),
        Bt_Bt1 = Bt + atan2(sinTheta, L - cosTheta);
    return [
      asin(x / sqrt(rho2)) * pi / Bt_Bt1,
      asin(1 - 2 * (theta - K2 * Bt - L * sinTheta + (1 + L2 - 2 * L * cosTheta) * Bt_Bt1) / B)
    ];
  };

  return forward;
}

/* harmony default export */ var hill = (function() {
  var K = 1,
      m = Object(d3_geo["geoProjectionMutator"])(hillRaw),
      p = m(K);

  p.ratio = function(_) {
    return arguments.length ? m(K = +_) : K;
  };

  return p
      .scale(167.774)
      .center([0, 18.67]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/sinuMollweide.js




var sinuMollweidePhi = 0.7109889596207567;

var sinuMollweideY = 0.0528035274542;

function sinuMollweideRaw(lambda, phi) {
  return phi > -sinuMollweidePhi
      ? (lambda = mollweideRaw(lambda, phi), lambda[1] += sinuMollweideY, lambda)
      : sinusoidalRaw(lambda, phi);
}

sinuMollweideRaw.invert = function(x, y) {
  return y > -sinuMollweidePhi
      ? mollweideRaw.invert(x, y - sinuMollweideY)
      : sinusoidalRaw.invert(x, y);
};

/* harmony default export */ var sinuMollweide = (function() {
  return Object(d3_geo["geoProjection"])(sinuMollweideRaw)
      .rotate([-20, -55])
      .scale(164.263)
      .center([0, -5.4036]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/homolosine.js






function homolosineRaw(lambda, phi) {
  return abs(phi) > sinuMollweidePhi
      ? (lambda = mollweideRaw(lambda, phi), lambda[1] -= phi > 0 ? sinuMollweideY : -sinuMollweideY, lambda)
      : sinusoidalRaw(lambda, phi);
}

homolosineRaw.invert = function(x, y) {
  return abs(y) > sinuMollweidePhi
      ? mollweideRaw.invert(x, y + (y > 0 ? sinuMollweideY : -sinuMollweideY))
      : sinusoidalRaw.invert(x, y);
};

/* harmony default export */ var homolosine = (function() {
  return Object(d3_geo["geoProjection"])(homolosineRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/index.js




function pointEqual(a, b) {
  return abs(a[0] - b[0]) < math_epsilon && abs(a[1] - b[1]) < math_epsilon;
}

function interpolateLine(coordinates, m) {
  var i = -1,
      n = coordinates.length,
      p0 = coordinates[0],
      p1,
      dx,
      dy,
      resampled = [];
  while (++i < n) {
    p1 = coordinates[i];
    dx = (p1[0] - p0[0]) / m;
    dy = (p1[1] - p0[1]) / m;
    for (var j = 0; j < m; ++j) resampled.push([p0[0] + j * dx, p0[1] + j * dy]);
    p0 = p1;
  }
  resampled.push(p1);
  return resampled;
}

function interpolateSphere(lobes) {
  var coordinates = [],
      lobe,
      lambda0, phi0, phi1,
      lambda2, phi2,
      i, n = lobes[0].length;

  // Northern Hemisphere
  for (i = 0; i < n; ++i) {
    lobe = lobes[0][i];
    lambda0 = lobe[0][0], phi0 = lobe[0][1], phi1 = lobe[1][1];
    lambda2 = lobe[2][0], phi2 = lobe[2][1];
    coordinates.push(interpolateLine([
      [lambda0 + math_epsilon, phi0 + math_epsilon],
      [lambda0 + math_epsilon, phi1 - math_epsilon],
      [lambda2 - math_epsilon, phi1 - math_epsilon],
      [lambda2 - math_epsilon, phi2 + math_epsilon]
    ], 30));
  }

  // Southern Hemisphere
  for (i = lobes[1].length - 1; i >= 0; --i) {
    lobe = lobes[1][i];
    lambda0 = lobe[0][0], phi0 = lobe[0][1], phi1 = lobe[1][1];
    lambda2 = lobe[2][0], phi2 = lobe[2][1];
    coordinates.push(interpolateLine([
      [lambda2 - math_epsilon, phi2 - math_epsilon],
      [lambda2 - math_epsilon, phi1 + math_epsilon],
      [lambda0 + math_epsilon, phi1 + math_epsilon],
      [lambda0 + math_epsilon, phi0 - math_epsilon]
    ], 30));
  }

  return {
    type: "Polygon",
    coordinates: [Object(src["a" /* merge */])(coordinates)]
  };
}

/* harmony default export */ var interrupted = (function(project, lobes) {
  var sphere = interpolateSphere(lobes);

  lobes = lobes.map(function(lobe) {
    return lobe.map(function(l) {
      return [
        [l[0][0] * radians, l[0][1] * radians],
        [l[1][0] * radians, l[1][1] * radians],
        [l[2][0] * radians, l[2][1] * radians]
      ];
    });
  });

  var bounds = lobes.map(function(lobe) {
    return lobe.map(function(l) {
      var x0 = project(l[0][0], l[0][1])[0],
          x1 = project(l[2][0], l[2][1])[0],
          y0 = project(l[1][0], l[0][1])[1],
          y1 = project(l[1][0], l[1][1])[1],
          t;
      if (y0 > y1) t = y0, y0 = y1, y1 = t;
      return [[x0, y0], [x1, y1]];
    });
  });

  function forward(lambda, phi) {
    var sign = phi < 0 ? -1 : +1, lobe = lobes[+(phi < 0)];
    for (var i = 0, n = lobe.length - 1; i < n && lambda > lobe[i][2][0]; ++i);
    var p = project(lambda - lobe[i][1][0], phi);
    p[0] += project(lobe[i][1][0], sign * phi > sign * lobe[i][0][1] ? lobe[i][0][1] : phi)[0];
    return p;
  }

  // Assumes mutually exclusive bounding boxes for lobes.
  if (project.invert) forward.invert = function(x, y) {
    var bound = bounds[+(y < 0)], lobe = lobes[+(y < 0)];
    for (var i = 0, n = bound.length; i < n; ++i) {
      var b = bound[i];
      if (b[0][0] <= x && x < b[1][0] && b[0][1] <= y && y < b[1][1]) {
        var p = project.invert(x - project(lobe[i][1][0], 0)[0], y);
        p[0] += lobe[i][1][0];
        return pointEqual(forward(p[0], p[1]), [x, y]) ? p : null;
      }
    }
  };

  var p = Object(d3_geo["geoProjection"])(forward),
      stream_ = p.stream;

  p.stream = function(stream) {
    var rotate = p.rotate(),
        rotateStream = stream_(stream),
        sphereStream = (p.rotate([0, 0]), stream_(stream));
    p.rotate(rotate);
    rotateStream.sphere = function() { Object(d3_geo["geoStream"])(sphere, sphereStream); };
    return rotateStream;
  };

  return p;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/boggs.js



var boggs_lobes = [[ // northern hemisphere
  [[-180,   0], [-100,  90], [ -40,   0]],
  [[ -40,   0], [  30,  90], [ 180,   0]]
], [ // southern hemisphere
  [[-180,   0], [-160, -90], [-100,   0]],
  [[-100,   0], [ -60, -90], [ -20,   0]],
  [[ -20,   0], [  20, -90], [  80,   0]],
  [[  80,   0], [ 140, -90], [ 180,   0]]
]];

/* harmony default export */ var interrupted_boggs = (function() {
  return interrupted(boggsRaw, boggs_lobes)
      .scale(160.857);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/homolosine.js



var homolosine_lobes = [[ // northern hemisphere
  [[-180,   0], [-100,  90], [ -40,   0]],
  [[ -40,   0], [  30,  90], [ 180,   0]]
], [ // southern hemisphere
  [[-180,   0], [-160, -90], [-100,   0]],
  [[-100,   0], [ -60, -90], [ -20,   0]],
  [[ -20,   0], [  20, -90], [  80,   0]],
  [[  80,   0], [ 140, -90], [ 180,   0]]
]];

/* harmony default export */ var interrupted_homolosine = (function() {
  return interrupted(homolosineRaw, homolosine_lobes)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/mollweide.js



var mollweide_lobes = [[ // northern hemisphere
  [[-180,   0], [-100,  90], [ -40,   0]],
  [[ -40,   0], [  30,  90], [ 180,   0]]
], [ // southern hemisphere
  [[-180,   0], [-160, -90], [-100,   0]],
  [[-100,   0], [ -60, -90], [ -20,   0]],
  [[ -20,   0], [  20, -90], [  80,   0]],
  [[  80,   0], [ 140, -90], [ 180,   0]]
]];

/* harmony default export */ var interrupted_mollweide = (function() {
  return interrupted(mollweideRaw, mollweide_lobes)
      .scale(169.529);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/mollweideHemispheres.js



var mollweideHemispheres_lobes = [[ // northern hemisphere
  [[-180,   0], [ -90,  90], [   0,   0]],
  [[   0,   0], [  90,  90], [ 180,   0]]
], [ // southern hemisphere
  [[-180,   0], [ -90, -90], [   0,   0]],
  [[   0,   0], [  90, -90], [ 180,   0]]
]];

/* harmony default export */ var mollweideHemispheres = (function() {
  return interrupted(mollweideRaw, mollweideHemispheres_lobes)
      .scale(169.529)
      .rotate([20, 0]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/sinuMollweide.js



var sinuMollweide_lobes = [[ // northern hemisphere
  [[-180,  35], [ -30,  90], [   0,  35]],
  [[   0,  35], [  30,  90], [ 180,  35]]
], [ // southern hemisphere
  [[-180, -10], [-102, -90], [ -65, -10]],
  [[ -65, -10], [   5, -90], [  77, -10]],
  [[  77, -10], [ 103, -90], [ 180, -10]]
]];

/* harmony default export */ var interrupted_sinuMollweide = (function() {
  return interrupted(sinuMollweideRaw, sinuMollweide_lobes)
      .rotate([-20, -55])
      .scale(164.263)
      .center([0, -5.4036]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/interrupted/sinusoidal.js



var sinusoidal_lobes = [[ // northern hemisphere
  [[-180,   0], [-110,  90], [ -40,   0]],
  [[ -40,   0], [   0,  90], [  40,   0]],
  [[  40,   0], [ 110,  90], [ 180,   0]]
], [ // southern hemisphere
  [[-180,   0], [-110, -90], [ -40,   0]],
  [[ -40,   0], [   0, -90], [  40,   0]],
  [[  40,   0], [ 110, -90], [ 180,   0]]
]];

/* harmony default export */ var interrupted_sinusoidal = (function() {
  return interrupted(sinusoidalRaw, sinusoidal_lobes)
      .scale(152.63)
      .rotate([-20, 0]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/kavrayskiy7.js



function kavrayskiy7Raw(lambda, phi) {
  return [3 / tau * lambda * sqrt(pi * pi / 3 - phi * phi), phi];
}

kavrayskiy7Raw.invert = function(x, y) {
  return [tau / 3 * x / sqrt(pi * pi / 3 - y * y), y];
};

/* harmony default export */ var kavrayskiy7 = (function() {
  return Object(d3_geo["geoProjection"])(kavrayskiy7Raw)
      .scale(158.837);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/lagrange.js



function lagrangeRaw(n) {

  function forward(lambda, phi) {
    if (abs(abs(phi) - halfPi) < math_epsilon) return [0, phi < 0 ? -2 : 2];
    var sinPhi = sin(phi),
        v = pow((1 + sinPhi) / (1 - sinPhi), n / 2),
        c = 0.5 * (v + 1 / v) + cos(lambda *= n);
    return [
      2 * sin(lambda) / c,
      (v - 1 / v) / c
    ];
  }

  forward.invert = function(x, y) {
    var y0 = abs(y);
    if (abs(y0 - 2) < math_epsilon) return x ? null : [0, sign(y) * halfPi];
    if (y0 > 2) return null;

    x /= 2, y /= 2;
    var x2 = x * x,
        y2 = y * y,
        t = 2 * y / (1 + x2 + y2); // tanh(nPhi)
    t = pow((1 + t) / (1 - t), 1 / n);
    return [
      atan2(2 * x, 1 - x2 - y2) / n,
      asin((t - 1) / (t + 1))
    ];
  };

  return forward;
}

/* harmony default export */ var lagrange = (function() {
  var n = 0.5,
      m = Object(d3_geo["geoProjectionMutator"])(lagrangeRaw),
      p = m(n);

  p.spacing = function(_) {
    return arguments.length ? m(n = +_) : n;
  };

  return p
      .scale(124.75);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/larrivee.js



var pi_sqrt2 = pi / sqrt2;

function larriveeRaw(lambda, phi) {
  return [
    lambda * (1 + sqrt(cos(phi))) / 2,
    phi / (cos(phi / 2) * cos(lambda / 6))
  ];
}

larriveeRaw.invert = function(x, y) {
  var x0 = abs(x),
      y0 = abs(y),
      lambda = math_epsilon,
      phi = halfPi;
  if (y0 < pi_sqrt2) phi *= y0 / pi_sqrt2;
  else lambda += 6 * acos(pi_sqrt2 / y0);
  for (var i = 0; i < 25; i++) {
    var sinPhi = sin(phi),
        sqrtcosPhi = sqrt(cos(phi)),
        sinPhi_2 = sin(phi / 2),
        cosPhi_2 = cos(phi / 2),
        sinLambda_6 = sin(lambda / 6),
        cosLambda_6 = cos(lambda / 6),
        f0 = 0.5 * lambda * (1 + sqrtcosPhi) - x0,
        f1 = phi / (cosPhi_2 * cosLambda_6) - y0,
        df0dPhi = sqrtcosPhi ? -0.25 * lambda * sinPhi / sqrtcosPhi : 0,
        df0dLambda = 0.5 * (1 + sqrtcosPhi),
        df1dPhi = (1 +0.5 * phi * sinPhi_2 / cosPhi_2) / (cosPhi_2 * cosLambda_6),
        df1dLambda = (phi / cosPhi_2) * (sinLambda_6 / 6) / (cosLambda_6 * cosLambda_6),
        denom = df0dPhi * df1dLambda - df1dPhi * df0dLambda,
        dPhi = (f0 * df1dLambda - f1 * df0dLambda) / denom,
        dLambda = (f1 * df0dPhi - f0 * df1dPhi) / denom;
    phi -= dPhi;
    lambda -= dLambda;
    if (abs(dPhi) < math_epsilon && abs(dLambda) < math_epsilon) break;
  }
  return [x < 0 ? -lambda : lambda, y < 0 ? -phi : phi];
};

/* harmony default export */ var larrivee = (function() {
  return Object(d3_geo["geoProjection"])(larriveeRaw)
      .scale(97.2672);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/laskowski.js



function laskowskiRaw(lambda, phi) {
  var lambda2 = lambda * lambda, phi2 = phi * phi;
  return [
    lambda * (0.975534 + phi2 * (-0.119161 + lambda2 * -0.0143059 + phi2 * -0.0547009)),
    phi * (1.00384 + lambda2 * (0.0802894 + phi2 * -0.02855 + lambda2 * 0.000199025) + phi2 * (0.0998909 + phi2 * -0.0491032))
  ];
}

laskowskiRaw.invert = function(x, y) {
  var lambda = sign(x) * pi,
      phi = y / 2,
      i = 50;
  do {
    var lambda2 = lambda * lambda,
        phi2 = phi * phi,
        lambdaPhi = lambda * phi,
        fx = lambda * (0.975534 + phi2 * (-0.119161 + lambda2 * -0.0143059 + phi2 * -0.0547009)) - x,
        fy = phi * (1.00384 + lambda2 * (0.0802894 + phi2 * -0.02855 + lambda2 * 0.000199025) + phi2 * (0.0998909 + phi2 * -0.0491032)) - y,
        deltaxDeltaLambda = 0.975534 - phi2 * (0.119161 + 3 * lambda2 * 0.0143059 + phi2 * 0.0547009),
        deltaxDeltaPhi = -lambdaPhi * (2 * 0.119161 + 4 * 0.0547009 * phi2 + 2 * 0.0143059 * lambda2),
        deltayDeltaLambda = lambdaPhi * (2 * 0.0802894 + 4 * 0.000199025 * lambda2 + 2 * -0.02855 * phi2),
        deltayDeltaPhi = 1.00384 + lambda2 * (0.0802894 + 0.000199025 * lambda2) + phi2 * (3 * (0.0998909 - 0.02855 * lambda2) - 5 * 0.0491032 * phi2),
        denominator = deltaxDeltaPhi * deltayDeltaLambda - deltayDeltaPhi * deltaxDeltaLambda,
        deltaLambda = (fy * deltaxDeltaPhi - fx * deltayDeltaPhi) / denominator,
        deltaPhi = (fx * deltayDeltaLambda - fy * deltaxDeltaLambda) / denominator;
    lambda -= deltaLambda, phi -= deltaPhi;
  } while ((abs(deltaLambda) > math_epsilon || abs(deltaPhi) > math_epsilon) && --i > 0);
  return i && [lambda, phi];
};

/* harmony default export */ var laskowski = (function() {
  return Object(d3_geo["geoProjection"])(laskowskiRaw)
      .scale(139.98);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/littrow.js



function littrowRaw(lambda, phi) {
  return [
    sin(lambda) / cos(phi),
    tan(phi) * cos(lambda)
  ];
}

littrowRaw.invert = function(x, y) {
  var x2 = x * x,
      y2 = y * y,
      y2_1 = y2 + 1,
      cosPhi = x
          ? sqrt1_2 * sqrt((y2_1 - sqrt(x2 * x2 + 2 * x2 * (y2 - 1) + y2_1 * y2_1)) / x2 + 1)
          : 1 / sqrt(y2_1);
  return [
    asin(x * cosPhi),
    sign(y) * acos(cosPhi)
  ];
};

/* harmony default export */ var littrow = (function() {
  return Object(d3_geo["geoProjection"])(littrowRaw)
      .scale(144.049)
      .clipAngle(90 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/loximuthal.js



function loximuthalRaw(phi0) {
  var cosPhi0 = cos(phi0),
      tanPhi0 = tan(quarterPi + phi0 / 2);

  function forward(lambda, phi) {
    var y = phi - phi0,
        x = abs(y) < math_epsilon ? lambda * cosPhi0
            : abs(x = quarterPi + phi / 2) < math_epsilon || abs(abs(x) - halfPi) < math_epsilon
            ? 0 : lambda * y / log(tan(x) / tanPhi0);
    return [x, y];
  }

  forward.invert = function(x, y) {
    var lambda,
        phi = y + phi0;
    return [
      abs(y) < math_epsilon ? x / cosPhi0
          : (abs(lambda = quarterPi + phi / 2) < math_epsilon || abs(abs(lambda) - halfPi) < math_epsilon) ? 0
          : x * log(tan(lambda) / tanPhi0) / y,
      phi
    ];
  };

  return forward;
}

/* harmony default export */ var loximuthal = (function() {
  return parallel1(loximuthalRaw)
      .parallel(40)
      .scale(158.837);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/miller.js



function millerRaw(lambda, phi) {
  return [lambda, 1.25 * log(tan(quarterPi + 0.4 * phi))];
}

millerRaw.invert = function(x, y) {
  return [x, 2.5 * atan(exp(0.8 * y)) - 0.625 * pi];
};

/* harmony default export */ var miller = (function() {
  return Object(d3_geo["geoProjection"])(millerRaw)
      .scale(108.318);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/modifiedStereographic.js



function modifiedStereographicRaw(C) {
  var m = C.length - 1;

  function forward(lambda, phi) {
    var cosPhi = cos(phi),
        k = 2 / (1 + cosPhi * cos(lambda)),
        zr = k * cosPhi * sin(lambda),
        zi = k * sin(phi),
        i = m,
        w = C[i],
        ar = w[0],
        ai = w[1],
        t;
    while (--i >= 0) {
      w = C[i];
      ar = w[0] + zr * (t = ar) - zi * ai;
      ai = w[1] + zr * ai + zi * t;
    }
    ar = zr * (t = ar) - zi * ai;
    ai = zr * ai + zi * t;
    return [ar, ai];
  }

  forward.invert = function(x, y) {
    var i = 20,
        zr = x,
        zi = y;
    do {
      var j = m,
          w = C[j],
          ar = w[0],
          ai = w[1],
          br = 0,
          bi = 0,
          t;

      while (--j >= 0) {
        w = C[j];
        br = ar + zr * (t = br) - zi * bi;
        bi = ai + zr * bi + zi * t;
        ar = w[0] + zr * (t = ar) - zi * ai;
        ai = w[1] + zr * ai + zi * t;
      }
      br = ar + zr * (t = br) - zi * bi;
      bi = ai + zr * bi + zi * t;
      ar = zr * (t = ar) - zi * ai - x;
      ai = zr * ai + zi * t - y;

      var denominator = br * br + bi * bi, deltar, deltai;
      zr -= deltar = (ar * br + ai * bi) / denominator;
      zi -= deltai = (ai * br - ar * bi) / denominator;
    } while (abs(deltar) + abs(deltai) > math_epsilon * math_epsilon && --i > 0);

    if (i) {
      var rho = sqrt(zr * zr + zi * zi),
          c = 2 * atan(rho * 0.5),
          sinc = sin(c);
      return [atan2(zr * sinc, rho * cos(c)), rho ? asin(zi * sinc / rho) : 0];
    }
  };

  return forward;
}

var alaska = [[0.9972523, 0], [0.0052513, -0.0041175], [0.0074606, 0.0048125], [-0.0153783, -0.1968253], [0.0636871, -0.1408027], [0.3660976, -0.2937382]],
    gs48 = [[0.98879, 0], [0, 0], [-0.050909, 0], [0, 0], [0.075528, 0]],
    gs50 = [[0.9842990, 0], [0.0211642, 0.0037608], [-0.1036018, -0.0575102], [-0.0329095, -0.0320119], [0.0499471, 0.1223335], [0.0260460, 0.0899805], [0.0007388, -0.1435792], [0.0075848, -0.1334108], [-0.0216473, 0.0776645], [-0.0225161, 0.0853673]],
    modifiedStereographic_miller = [[0.9245, 0], [0, 0], [0.01943, 0]],
    lee = [[0.721316, 0], [0, 0], [-0.00881625, -0.00617325]];

function modifiedStereographicAlaska() {
  return modifiedStereographic(alaska, [152, -64])
      .scale(1500)
      .center([-160.908, 62.4864])
      .clipAngle(25);
}

function modifiedStereographicGs48() {
  return modifiedStereographic(gs48, [95, -38])
      .scale(1000)
      .clipAngle(55)
      .center([-96.5563, 38.8675]);
}

function modifiedStereographicGs50() {
  return modifiedStereographic(gs50, [120, -45])
      .scale(359.513)
      .clipAngle(55)
      .center([-117.474, 53.0628]);
}

function modifiedStereographicMiller() {
  return modifiedStereographic(modifiedStereographic_miller, [-20, -18])
      .scale(209.091)
      .center([20, 16.7214])
      .clipAngle(82);
}

function modifiedStereographicLee() {
  return modifiedStereographic(lee, [165, 10])
      .scale(250)
      .clipAngle(130)
      .center([-165, -10]);
}

function modifiedStereographic(coefficients, rotate) {
  var p = Object(d3_geo["geoProjection"])(modifiedStereographicRaw(coefficients)).rotate(rotate).clipAngle(90),
      r = Object(d3_geo["geoRotation"])(rotate),
      center = p.center;

  delete p.rotate;

  p.center = function(_) {
    return arguments.length ? center(r(_)) : r.invert(center());
  };

  return p;
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/mtFlatPolarParabolic.js



var sqrt6 = sqrt(6),
    sqrt7 = sqrt(7);

function mtFlatPolarParabolicRaw(lambda, phi) {
  var theta = asin(7 * sin(phi) / (3 * sqrt6));
  return [
    sqrt6 * lambda * (2 * cos(2 * theta / 3) - 1) / sqrt7,
    9 * sin(theta / 3) / sqrt7
  ];
}

mtFlatPolarParabolicRaw.invert = function(x, y) {
  var theta = 3 * asin(y * sqrt7 / 9);
  return [
    x * sqrt7 / (sqrt6 * (2 * cos(2 * theta / 3) - 1)),
    asin(sin(theta) * 3 * sqrt6 / 7)
  ];
};

/* harmony default export */ var mtFlatPolarParabolic = (function() {
  return Object(d3_geo["geoProjection"])(mtFlatPolarParabolicRaw)
      .scale(164.859);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/mtFlatPolarQuartic.js



function mtFlatPolarQuarticRaw(lambda, phi) {
  var k = (1 + sqrt1_2) * sin(phi),
      theta = phi;
  for (var i = 0, delta; i < 25; i++) {
    theta -= delta = (sin(theta / 2) + sin(theta) - k) / (0.5 * cos(theta / 2) + cos(theta));
    if (abs(delta) < math_epsilon) break;
  }
  return [
    lambda * (1 + 2 * cos(theta) / cos(theta / 2)) / (3 * sqrt2),
    2 * sqrt(3) * sin(theta / 2) / sqrt(2 + sqrt2)
  ];
}

mtFlatPolarQuarticRaw.invert = function(x, y) {
  var sinTheta_2 = y * sqrt(2 + sqrt2) / (2 * sqrt(3)),
      theta = 2 * asin(sinTheta_2);
  return [
    3 * sqrt2 * x / (1 + 2 * cos(theta) / cos(theta / 2)),
    asin((sinTheta_2 + sin(theta)) / (1 + sqrt1_2))
  ];
};

/* harmony default export */ var mtFlatPolarQuartic = (function() {
  return Object(d3_geo["geoProjection"])(mtFlatPolarQuarticRaw)
      .scale(188.209);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/mtFlatPolarSinusoidal.js



function mtFlatPolarSinusoidalRaw(lambda, phi) {
  var A = sqrt(6 / (4 + pi)),
      k = (1 + pi / 4) * sin(phi),
      theta = phi / 2;
  for (var i = 0, delta; i < 25; i++) {
    theta -= delta = (theta / 2 + sin(theta) - k) / (0.5 + cos(theta));
    if (abs(delta) < math_epsilon) break;
  }
  return [
    A * (0.5 + cos(theta)) * lambda / 1.5,
    A * theta
  ];
}

mtFlatPolarSinusoidalRaw.invert = function(x, y) {
  var A = sqrt(6 / (4 + pi)),
      theta = y / A;
  if (abs(abs(theta) - halfPi) < math_epsilon) theta = theta < 0 ? -halfPi : halfPi;
  return [
    1.5 * x / (A * (0.5 + cos(theta))),
    asin((theta / 2 + sin(theta)) / (1 + pi / 4))
  ];
};

/* harmony default export */ var mtFlatPolarSinusoidal = (function() {
  return Object(d3_geo["geoProjection"])(mtFlatPolarSinusoidalRaw)
      .scale(166.518);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/naturalEarth.js



function naturalEarthRaw(lambda, phi) {
  var phi2 = phi * phi, phi4 = phi2 * phi2;
  return [
    lambda * (0.8707 - 0.131979 * phi2 + phi4 * (-0.013791 + phi4 * (0.003971 * phi2 - 0.001529 * phi4))),
    phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4)))
  ];
}

naturalEarthRaw.invert = function(x, y) {
  var phi = y, i = 25, delta;
  do {
    var phi2 = phi * phi, phi4 = phi2 * phi2;
    phi -= delta = (phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4))) - y) /
        (1.007226 + phi2 * (0.015085 * 3 + phi4 * (-0.044475 * 7 + 0.028874 * 9 * phi2 - 0.005916 * 11 * phi4)));
  } while (abs(delta) > math_epsilon && --i > 0);
  return [
    x / (0.8707 + (phi2 = phi * phi) * (-0.131979 + phi2 * (-0.013791 + phi2 * phi2 * phi2 * (0.003971 - 0.001529 * phi2)))),
    phi
  ];
};

/* harmony default export */ var naturalEarth = (function() {
  return Object(d3_geo["geoProjection"])(naturalEarthRaw)
      .scale(175.295);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/nellHammer.js



function nellHammerRaw(lambda, phi) {
  return [
    lambda * (1 + cos(phi)) / 2,
    2 * (phi - tan(phi / 2))
  ];
}

nellHammerRaw.invert = function(x, y) {
  var p = y / 2;
  for (var i = 0, delta = Infinity; i < 10 && abs(delta) > math_epsilon; ++i) {
    var c = cos(y / 2);
    y -= delta = (y - tan(y / 2) - p) / (1 - 0.5 / (c * c));
  }
  return [
    2 * x / (1 + cos(y)),
    y
  ];
};

/* harmony default export */ var nellHammer = (function() {
  return Object(d3_geo["geoProjection"])(nellHammerRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/patterson.js



// Based on Java implementation by Bojan Savric.
// https://github.com/OSUCartography/JMapProjLib/blob/master/src/com/jhlabs/map/proj/PattersonProjection.java

var pattersonK1 = 1.0148,
    pattersonK2 = 0.23185,
    pattersonK3 = -0.14499,
    pattersonK4 = 0.02406,
    pattersonC1 = pattersonK1,
    pattersonC2 = 5 * pattersonK2,
    pattersonC3 = 7 * pattersonK3,
    pattersonC4 = 9 * pattersonK4,
    pattersonYmax = 1.790857183;

function pattersonRaw(lambda, phi) {
  var phi2 = phi * phi;
  return [
    lambda,
    phi * (pattersonK1 + phi2 * phi2 * (pattersonK2 + phi2 * (pattersonK3 + pattersonK4 * phi2)))
  ];
}

pattersonRaw.invert = function(x, y) {
  if (y > pattersonYmax) y = pattersonYmax;
  else if (y < -pattersonYmax) y = -pattersonYmax;
  var yc = y, delta;

  do { // Newton-Raphson
    var y2 = yc * yc;
    yc -= delta = ((yc * (pattersonK1 + y2 * y2 * (pattersonK2 + y2 * (pattersonK3 + pattersonK4 * y2)))) - y) / (pattersonC1 + y2 * y2 * (pattersonC2 + y2 * (pattersonC3 + pattersonC4 * y2)));
  } while (abs(delta) > math_epsilon);

  return [x, yc];
};

/* harmony default export */ var patterson = (function() {
  return Object(d3_geo["geoProjection"])(pattersonRaw)
      .scale(139.319);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyconic.js



function polyconicRaw(lambda, phi) {
  if (abs(phi) < math_epsilon) return [lambda, 0];
  var tanPhi = tan(phi),
      k = lambda * sin(phi);
  return [
    sin(k) / tanPhi,
    phi + (1 - cos(k)) / tanPhi
  ];
}

polyconicRaw.invert = function(x, y) {
  if (abs(y) < math_epsilon) return [x, 0];
  var k = x * x + y * y,
      phi = y * 0.5,
      i = 10, delta;
  do {
    var tanPhi = tan(phi),
        secPhi = 1 / cos(phi),
        j = k - 2 * y * phi + phi * phi;
    phi -= delta = (tanPhi * j + 2 * (phi - y)) / (2 + j * secPhi * secPhi + 2 * (phi - y) * tanPhi);
  } while (abs(delta) > math_epsilon && --i > 0);
  tanPhi = tan(phi);
  return [
    (abs(y) < abs(phi + 1 / tanPhi) ? asin(x * tanPhi) : sign(x) * (acos(abs(x * tanPhi)) + halfPi)) / sin(phi),
    phi
  ];
};

/* harmony default export */ var polyconic = (function() {
  return Object(d3_geo["geoProjection"])(polyconicRaw)
      .scale(103.74);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyhedral/matrix.js


// Note: 6-element arrays are used to denote the 3x3 affine transform matrix:
// [a, b, c,
//  d, e, f,
//  0, 0, 1] - this redundant row is left out.

// Transform matrix for [a0, a1] -> [b0, b1].
/* harmony default export */ var matrix = (function(a, b) {
  var u = subtract(a[1], a[0]),
      v = subtract(b[1], b[0]),
      phi = matrix_angle(u, v),
      s = matrix_length(u) / matrix_length(v);

  return multiply([
    1, 0, a[0][0],
    0, 1, a[0][1]
  ], multiply([
    s, 0, 0,
    0, s, 0
  ], multiply([
    cos(phi), sin(phi), 0,
    -sin(phi), cos(phi), 0
  ], [
    1, 0, -b[0][0],
    0, 1, -b[0][1]
  ])));
});

// Inverts a transform matrix.
function inverse(m) {
  var k = 1 / (m[0] * m[4] - m[1] * m[3]);
  return [
    k * m[4], -k * m[1], k * (m[1] * m[5] - m[2] * m[4]),
    -k * m[3], k * m[0], k * (m[2] * m[3] - m[0] * m[5])
  ];
}

// Multiplies two 3x2 matrices.
function multiply(a, b) {
  return [
    a[0] * b[0] + a[1] * b[3],
    a[0] * b[1] + a[1] * b[4],
    a[0] * b[2] + a[1] * b[5] + a[2],
    a[3] * b[0] + a[4] * b[3],
    a[3] * b[1] + a[4] * b[4],
    a[3] * b[2] + a[4] * b[5] + a[5]
  ];
}

// Subtracts 2D vectors.
function subtract(a, b) {
  return [a[0] - b[0], a[1] - b[1]];
}

// Magnitude of a 2D vector.
function matrix_length(v) {
  return sqrt(v[0] * v[0] + v[1] * v[1]);
}

// Angle between two 2D vectors.
function matrix_angle(a, b) {
  return atan2(a[0] * b[1] - a[1] * b[0], a[0] * b[0] + a[1] * b[1]);
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyhedral/index.js




// Creates a polyhedral projection.
//  * root: a spanning tree of polygon faces.  Nodes are automatically
//    augmented with a transform matrix.
//  * face: a function that returns the appropriate node for a given {lambda, phi}
//    point (radians).
//  * r: rotation angle for final polyhedral net.  Defaults to -pi / 6 (for
//    butterflies).
/* harmony default export */ var polyhedral = (function(root, face, r) {

  r = r == null ? -pi / 6 : r; // TODO automate

  recurse(root, {transform: [
    cos(r), sin(r), 0,
    -sin(r), cos(r), 0
  ]});

  function recurse(node, parent) {
    node.edges = faceEdges(node.face);
    // Find shared edge.
    if (parent.face) {
      var shared = node.shared = sharedEdge(node.face, parent.face),
          m = matrix(shared.map(parent.project), shared.map(node.project));
      node.transform = parent.transform ? multiply(parent.transform, m) : m;
      // Replace shared edge in parent edges array.
      var edges = parent.edges;
      for (var i = 0, n = edges.length; i < n; ++i) {
        if (polyhedral_pointEqual(shared[0], edges[i][1]) && polyhedral_pointEqual(shared[1], edges[i][0])) edges[i] = node;
        if (polyhedral_pointEqual(shared[0], edges[i][0]) && polyhedral_pointEqual(shared[1], edges[i][1])) edges[i] = node;
      }
      edges = node.edges;
      for (i = 0, n = edges.length; i < n; ++i) {
        if (polyhedral_pointEqual(shared[0], edges[i][0]) && polyhedral_pointEqual(shared[1], edges[i][1])) edges[i] = parent;
        if (polyhedral_pointEqual(shared[0], edges[i][1]) && polyhedral_pointEqual(shared[1], edges[i][0])) edges[i] = parent;
      }
    } else {
      node.transform = parent.transform;
    }
    if (node.children) {
      node.children.forEach(function(child) {
        recurse(child, node);
      });
    }
    return node;
  }

  function forward(lambda, phi) {
    var node = face(lambda, phi),
        point = node.project([lambda * degrees, phi * degrees]),
        t;
    if (t = node.transform) {
      return [
        t[0] * point[0] + t[1] * point[1] + t[2],
        -(t[3] * point[0] + t[4] * point[1] + t[5])
      ];
    }
    point[1] = -point[1];
    return point;
  }

  // Naive inverse!  A faster solution would use bounding boxes, or even a
  // polygonal quadtree.
  if (hasInverse(root)) forward.invert = function(x, y) {
    var coordinates = faceInvert(root, [x, -y]);
    return coordinates && (coordinates[0] *= radians, coordinates[1] *= radians, coordinates);
  };

  function faceInvert(node, coordinates) {
    var invert = node.project.invert,
        t = node.transform,
        point = coordinates;
    if (t) {
      t = inverse(t);
      point = [
        t[0] * point[0] + t[1] * point[1] + t[2],
        (t[3] * point[0] + t[4] * point[1] + t[5])
      ];
    }
    if (invert && node === faceDegrees(p = invert(point))) return p;
    var p,
        children = node.children;
    for (var i = 0, n = children && children.length; i < n; ++i) {
      if (p = faceInvert(children[i], coordinates)) return p;
    }
  }

  function faceDegrees(coordinates) {
    return face(coordinates[0] * radians, coordinates[1] * radians);
  }

  var proj = Object(d3_geo["geoProjection"])(forward),
      stream_ = proj.stream;

  proj.stream = function(stream) {
    var rotate = proj.rotate(),
        rotateStream = stream_(stream),
        sphereStream = (proj.rotate([0, 0]), stream_(stream));
    proj.rotate(rotate);
    rotateStream.sphere = function() {
      sphereStream.polygonStart();
      sphereStream.lineStart();
      outline(sphereStream, root);
      sphereStream.lineEnd();
      sphereStream.polygonEnd();
    };
    return rotateStream;
  };

  return proj;
});

function outline(stream, node, parent) {
  var point,
      edges = node.edges,
      n = edges.length,
      edge,
      multiPoint = {type: "MultiPoint", coordinates: node.face},
      notPoles = node.face.filter(function(d) { return abs(d[1]) !== 90; }),
      b = Object(d3_geo["geoBounds"])({type: "MultiPoint", coordinates: notPoles}),
      inside = false,
      j = -1,
      dx = b[1][0] - b[0][0];
  // TODO
  var c = dx === 180 || dx === 360
      ? [(b[0][0] + b[1][0]) / 2, (b[0][1] + b[1][1]) / 2]
      : Object(d3_geo["geoCentroid"])(multiPoint);
  // First find the shared edge…
  if (parent) while (++j < n) {
    if (edges[j] === parent) break;
  }
  ++j;
  for (var i = 0; i < n; ++i) {
    edge = edges[(i + j) % n];
    if (Array.isArray(edge)) {
      if (!inside) {
        stream.point((point = Object(d3_geo["geoInterpolate"])(edge[0], c)(math_epsilon))[0], point[1]);
        inside = true;
      }
      stream.point((point = Object(d3_geo["geoInterpolate"])(edge[1], c)(math_epsilon))[0], point[1]);
    } else {
      inside = false;
      if (edge !== parent) outline(stream, edge, node);
    }
  }
}

// Tests equality of two spherical points.
function polyhedral_pointEqual(a, b) {
  return a && b && a[0] === b[0] && a[1] === b[1];
}

// Finds a shared edge given two clockwise polygons.
function sharedEdge(a, b) {
  var x, y, n = a.length, found = null;
  for (var i = 0; i < n; ++i) {
    x = a[i];
    for (var j = b.length; --j >= 0;) {
      y = b[j];
      if (x[0] === y[0] && x[1] === y[1]) {
        if (found) return [found, x];
        found = x;
      }
    }
  }
}

// Converts an array of n face vertices to an array of n + 1 edges.
function faceEdges(face) {
  var n = face.length,
      edges = [];
  for (var a = face[n - 1], i = 0; i < n; ++i) edges.push([a, a = face[i]]);
  return edges;
}

function hasInverse(node) {
  return node.project.invert || node.children && node.children.some(hasInverse);
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyhedral/octahedron.js
// TODO generate on-the-fly to avoid external modification.
var octahedron = [
  [0, 90],
  [-90, 0], [0, 0], [90, 0], [180, 0],
  [0, -90]
];

/* harmony default export */ var polyhedral_octahedron = ([
  [0, 2, 1],
  [0, 3, 2],
  [5, 1, 2],
  [5, 2, 3],
  [0, 1, 4],
  [0, 4, 3],
  [5, 4, 1],
  [5, 3, 4]
].map(function(face) {
  return face.map(function(i) {
    return octahedron[i];
  });
}));

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyhedral/butterfly.js





/* harmony default export */ var butterfly = (function(faceProjection) {

  faceProjection = faceProjection || function(face) {
    var c = Object(d3_geo["geoCentroid"])({type: "MultiPoint", coordinates: face});
    return Object(d3_geo["geoGnomonic"])().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]);
  };

  var faces = polyhedral_octahedron.map(function(face) {
    return {face: face, project: faceProjection(face)};
  });

  [-1, 0, 0, 1, 0, 1, 4, 5].forEach(function(d, i) {
    var node = faces[d];
    node && (node.children || (node.children = [])).push(faces[i]);
  });

  return polyhedral(faces[0], function(lambda, phi) {
        return faces[lambda < -pi / 2 ? phi < 0 ? 6 : 4
            : lambda < 0 ? phi < 0 ? 2 : 0
            : lambda < pi / 2 ? phi < 0 ? 3 : 1
            : phi < 0 ? 7 : 5];
      })
      .scale(101.858)
      .center([0, 45]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyhedral/collignon.js






var kx = 2 / sqrt(3);

function collignonK(a, b) {
  var p = collignonRaw(a, b);
  return [p[0] * kx, p[1]];
}

collignonK.invert = function(x,y) {
  return collignonRaw.invert(x / kx, y);
};

/* harmony default export */ var polyhedral_collignon = (function(faceProjection) {

  faceProjection = faceProjection || function(face) {
    var c = Object(d3_geo["geoCentroid"])({type: "MultiPoint", coordinates: face});
    return Object(d3_geo["geoProjection"])(collignonK).translate([0, 0]).scale(1).rotate(c[1] > 0 ? [-c[0], 0] : [180 - c[0], 180]);
  };

  var faces = polyhedral_octahedron.map(function(face) {
    return {face: face, project: faceProjection(face)};
  });

  [-1, 0, 0, 1, 0, 1, 4, 5].forEach(function(d, i) {
    var node = faces[d];
    node && (node.children || (node.children = [])).push(faces[i]);
  });

  return polyhedral(faces[0], function(lambda, phi) {
        return faces[lambda < -pi / 2 ? phi < 0 ? 6 : 4
            : lambda < 0 ? phi < 0 ? 2 : 0
            : lambda < pi / 2 ? phi < 0 ? 3 : 1
            : phi < 0 ? 7 : 5];
      })
      .scale(121.906)
      .center([0, 48.5904]);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/polyhedral/waterman.js





/* harmony default export */ var waterman = (function(faceProjection) {

  faceProjection = faceProjection || function(face) {
    var c = face.length === 6 ? Object(d3_geo["geoCentroid"])({type: "MultiPoint", coordinates: face}) : face[0];
    return Object(d3_geo["geoGnomonic"])().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]);
  };

  var w5 = polyhedral_octahedron.map(function(face) {
    var xyz = face.map(waterman_cartesian),
        n = xyz.length,
        a = xyz[n - 1],
        b,
        hexagon = [];
    for (var i = 0; i < n; ++i) {
      b = xyz[i];
      hexagon.push(spherical([
        a[0] * 0.9486832980505138 + b[0] * 0.31622776601683794,
        a[1] * 0.9486832980505138 + b[1] * 0.31622776601683794,
        a[2] * 0.9486832980505138 + b[2] * 0.31622776601683794
      ]), spherical([
        b[0] * 0.9486832980505138 + a[0] * 0.31622776601683794,
        b[1] * 0.9486832980505138 + a[1] * 0.31622776601683794,
        b[2] * 0.9486832980505138 + a[2] * 0.31622776601683794
      ]));
      a = b;
    }
    return hexagon;
  });

  var cornerNormals = [];

  var parents = [-1, 0, 0, 1, 0, 1, 4, 5];

  w5.forEach(function(hexagon, j) {
    var face = polyhedral_octahedron[j],
        n = face.length,
        normals = cornerNormals[j] = [];
    for (var i = 0; i < n; ++i) {
      w5.push([
        face[i],
        hexagon[(i * 2 + 2) % (2 * n)],
        hexagon[(i * 2 + 1) % (2 * n)]
      ]);
      parents.push(j);
      normals.push(cross(
        waterman_cartesian(hexagon[(i * 2 + 2) % (2 * n)]),
        waterman_cartesian(hexagon[(i * 2 + 1) % (2 * n)])
      ));
    }
  });

  var faces = w5.map(function(face) {
    return {
      project: faceProjection(face),
      face: face
    };
  });

  parents.forEach(function(d, i) {
    var parent = faces[d];
    parent && (parent.children || (parent.children = [])).push(faces[i]);
  });

  function face(lambda, phi) {
    var cosphi = cos(phi),
        p = [cosphi * cos(lambda), cosphi * sin(lambda), sin(phi)];

    var hexagon = lambda < -pi / 2 ? phi < 0 ? 6 : 4
        : lambda < 0 ? phi < 0 ? 2 : 0
        : lambda < pi / 2 ? phi < 0 ? 3 : 1
        : phi < 0 ? 7 : 5;

    var n = cornerNormals[hexagon];

    return faces[dot(n[0], p) < 0 ? 8 + 3 * hexagon
        : dot(n[1], p) < 0 ? 8 + 3 * hexagon + 1
        : dot(n[2], p) < 0 ? 8 + 3 * hexagon + 2
        : hexagon];
  }

  return polyhedral(faces[0], face)
      .scale(110.625)
      .center([0,45]);
});

function dot(a, b) {
  for (var i = 0, n = a.length, s = 0; i < n; ++i) s += a[i] * b[i];
  return s;
}

function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0]
  ];
}

// Converts 3D Cartesian to spherical coordinates (degrees).
function spherical(cartesian) {
  return [
    atan2(cartesian[1], cartesian[0]) * degrees,
    asin(max(-1, min(1, cartesian[2]))) * degrees
  ];
}

// Converts spherical coordinates (degrees) to 3D Cartesian.
function waterman_cartesian(coordinates) {
  var lambda = coordinates[0] * radians,
      phi = coordinates[1] * radians,
      cosphi = cos(phi);
  return [
    cosphi * cos(lambda),
    cosphi * sin(lambda),
    sin(phi)
  ];
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/noop.js
/* harmony default export */ var noop = (function() {});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/project/clockwise.js
/* harmony default export */ var clockwise = (function(ring) {
  if ((n = ring.length) < 4) return false;
  var i = 0,
      n,
      area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
  while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
  return area <= 0;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/project/contains.js
/* harmony default export */ var contains = (function(ring, point) {
  var x = point[0],
      y = point[1],
      contains = false;
  for (var i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
    var pi = ring[i], xi = pi[0], yi = pi[1],
        pj = ring[j], xj = pj[0], yj = pj[1];
    if (((yi > y) ^ (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) contains = !contains;
  }
  return contains;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/project/index.js





/* harmony default export */ var src_project = (function(object, projection) {
  var stream = projection.stream, project;
  if (!stream) throw new Error("invalid projection");
  switch (object && object.type) {
    case "Feature": project = projectFeature; break;
    case "FeatureCollection": project = projectFeatureCollection; break;
    default: project = projectGeometry; break;
  }
  return project(object, stream);
});

function projectFeatureCollection(o, stream) {
  return {
    type: "FeatureCollection",
    features: o.features.map(function(f) {
      return projectFeature(f, stream);
    })
  };
}

function projectFeature(o, stream) {
  return {
    type: "Feature",
    id: o.id,
    properties: o.properties,
    geometry: projectGeometry(o.geometry, stream)
  };
}

function projectGeometryCollection(o, stream) {
  return {
    type: "GeometryCollection",
    geometries: o.geometries.map(function(o) {
      return projectGeometry(o, stream);
    })
  };
}

function projectGeometry(o, stream) {
  if (!o) return null;
  if (o.type === "GeometryCollection") return projectGeometryCollection(o, stream);
  var sink;
  switch (o.type) {
    case "Point": sink = sinkPoint; break;
    case "MultiPoint": sink = sinkPoint; break;
    case "LineString": sink = sinkLine; break;
    case "MultiLineString": sink = sinkLine; break;
    case "Polygon": sink = sinkPolygon; break;
    case "MultiPolygon": sink = sinkPolygon; break;
    case "Sphere": sink = sinkPolygon; break;
    default: return null;
  }
  Object(d3_geo["geoStream"])(o, stream(sink));
  return sink.result();
}

var project_points = [],
    lines = [];

var sinkPoint = {
  point: function(x, y) {
    project_points.push([x, y]);
  },
  result: function() {
    var result = !project_points.length ? null
        : project_points.length < 2 ? {type: "Point", coordinates: project_points[0]}
        : {type: "MultiPoint", coordinates: project_points};
    project_points = [];
    return result;
  }
};

var sinkLine = {
  lineStart: noop,
  point: function(x, y) {
    project_points.push([x, y]);
  },
  lineEnd: function() {
    if (project_points.length) lines.push(project_points), project_points = [];
  },
  result: function() {
    var result = !lines.length ? null
        : lines.length < 2 ? {type: "LineString", coordinates: lines[0]}
        : {type: "MultiLineString", coordinates: lines};
    lines = [];
    return result;
  }
};

var sinkPolygon = {
  polygonStart: noop,
  lineStart: noop,
  point: function(x, y) {
    project_points.push([x, y]);
  },
  lineEnd: function() {
    var n = project_points.length;
    if (n) {
      do project_points.push(project_points[0].slice()); while (++n < 4);
      lines.push(project_points), project_points = [];
    }
  },
  polygonEnd: noop,
  result: function() {
    if (!lines.length) return null;
    var polygons = [],
        holes = [];

    // https://github.com/d3/d3/issues/1558
    lines.forEach(function(ring) {
      if (clockwise(ring)) polygons.push([ring]);
      else holes.push(ring);
    });

    holes.forEach(function(hole) {
      var point = hole[0];
      polygons.some(function(polygon) {
        if (contains(polygon[0], point)) {
          polygon.push(hole);
          return true;
        }
      }) || polygons.push([hole]);
    });

    lines = [];

    return !polygons.length ? null
        : polygons.length > 1 ? {type: "MultiPolygon", coordinates: polygons}
        : {type: "Polygon", coordinates: polygons[0]};
  }
};

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/quincuncial/index.js



/* harmony default export */ var quincuncial = (function(project) {
  var dx = project(halfPi, 0)[0] - project(-halfPi, 0)[0];

  function projectQuincuncial(lambda, phi) {
    var t = abs(lambda) < halfPi,
        p = project(t ? lambda : lambda > 0 ? lambda - pi : lambda + pi, phi),
        x = (p[0] - p[1]) * sqrt1_2,
        y = (p[0] + p[1]) * sqrt1_2;
    if (t) return [x, y];
    var d = dx * sqrt1_2,
        s = x > 0 ^ y > 0 ? -1 : 1;
    return [s * x - sign(y) * d, s * y - sign(x) * d];
  }

  if (project.invert) projectQuincuncial.invert = function(x0, y0) {
    var x = (x0 + y0) * sqrt1_2,
        y = (y0 - x0) * sqrt1_2,
        t = abs(x) < 0.5 * dx && abs(y) < 0.5 * dx;

    if (!t) {
      var d = dx * sqrt1_2,
          s = x > 0 ^ y > 0 ? -1 : 1,
          x1 = -s * (x0 + (y > 0 ? 1 : -1) * d),
          y1 = -s * (y0 + (x > 0 ? 1 : -1) * d);
      x = (-x1 - y1) * sqrt1_2;
      y = (x1 - y1) * sqrt1_2;
    }

    var p = project.invert(x, y);
    if (!t) p[0] += x > 0 ? pi : -pi;
    return p;
  };

  return Object(d3_geo["geoProjection"])(projectQuincuncial)
      .rotate([-90, -90, 45])
      .clipAngle(180 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/quincuncial/gringorten.js



/* harmony default export */ var quincuncial_gringorten = (function() {
  return quincuncial(gringortenRaw)
      .scale(176.423);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/quincuncial/peirce.js



/* harmony default export */ var peirce = (function() {
  return quincuncial(guyouRaw)
      .scale(111.48);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/quantize.js
/* harmony default export */ var quantize = (function(o, digits) {
  if (!(0 <= (digits = +digits) && digits <= 20)) throw new Error("invalid digits");

  function quantizePoint(coordinates) {
    coordinates[0] = +coordinates[0].toFixed(digits);
    coordinates[1] = +coordinates[1].toFixed(digits);
  }

  function quantizePoints(coordinates) {
    coordinates.forEach(quantizePoint);
  }

  function quantizePolygon(coordinates) {
    coordinates.forEach(quantizePoints);
  }

  function quantizeGeometry(o) {
    if (o) switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(quantizeGeometry); break;
      case "Point": quantizePoint(o.coordinates); break;
      case "MultiPoint": case "LineString": quantizePoints(o.coordinates); break;
      case "MultiLineString": case "Polygon": quantizePolygon(o.coordinates); break;
      case "MultiPolygon": o.coordinates.forEach(quantizePolygon); break;
      default: return;
    }
  }

  function quantizeFeature(o) {
    quantizeGeometry(o.geometry);
  }

  if (o) switch (o.type) {
    case "Feature": quantizeFeature(o); break;
    case "FeatureCollection": o.features.forEach(quantizeFeature); break;
    default: quantizeGeometry(o); break;
  }

  return o;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/rectangularPolyconic.js



function rectangularPolyconicRaw(phi0) {
  var sinPhi0 = sin(phi0);

  function forward(lambda, phi) {
    var A = sinPhi0 ? tan(lambda * sinPhi0 / 2) / sinPhi0 : lambda / 2;
    if (!phi) return [2 * A, -phi0];
    var E = 2 * atan(A * sin(phi)),
        cotPhi = 1 / tan(phi);
    return [
      sin(E) * cotPhi,
      phi + (1 - cos(E)) * cotPhi - phi0
    ];
  }

  // TODO return null for points outside outline.
  forward.invert = function(x, y) {
    if (abs(y += phi0) < math_epsilon) return [sinPhi0 ? 2 * atan(sinPhi0 * x / 2) / sinPhi0 : x, 0];
    var k = x * x + y * y,
        phi = 0,
        i = 10, delta;
    do {
      var tanPhi = tan(phi),
          secPhi = 1 / cos(phi),
          j = k - 2 * y * phi + phi * phi;
      phi -= delta = (tanPhi * j + 2 * (phi - y)) / (2 + j * secPhi * secPhi + 2 * (phi - y) * tanPhi);
    } while (abs(delta) > math_epsilon && --i > 0);
    var E = x * (tanPhi = tan(phi)),
        A = tan(abs(y) < abs(phi + 1 / tanPhi) ? asin(E) * 0.5 : acos(E) * 0.5 + pi / 4) / sin(phi);
    return [
      sinPhi0 ? 2 * atan(sinPhi0 * A) / sinPhi0 : 2 * A,
      phi
    ];
  };

  return forward;
}

/* harmony default export */ var rectangularPolyconic = (function() {
  return parallel1(rectangularPolyconicRaw)
      .scale(131.215);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/robinson.js



var robinson_K = [
  [0.9986, -0.062],
  [1.0000, 0.0000],
  [0.9986, 0.0620],
  [0.9954, 0.1240],
  [0.9900, 0.1860],
  [0.9822, 0.2480],
  [0.9730, 0.3100],
  [0.9600, 0.3720],
  [0.9427, 0.4340],
  [0.9216, 0.4958],
  [0.8962, 0.5571],
  [0.8679, 0.6176],
  [0.8350, 0.6769],
  [0.7986, 0.7346],
  [0.7597, 0.7903],
  [0.7186, 0.8435],
  [0.6732, 0.8936],
  [0.6213, 0.9394],
  [0.5722, 0.9761],
  [0.5322, 1.0000]
];

robinson_K.forEach(function(d) {
  d[1] *= 1.0144;
});

function robinsonRaw(lambda, phi) {
  var i = min(18, abs(phi) * 36 / pi),
      i0 = floor(i),
      di = i - i0,
      ax = (k = robinson_K[i0])[0],
      ay = k[1],
      bx = (k = robinson_K[++i0])[0],
      by = k[1],
      cx = (k = robinson_K[min(19, ++i0)])[0],
      cy = k[1],
      k;
  return [
    lambda * (bx + di * (cx - ax) / 2 + di * di * (cx - 2 * bx + ax) / 2),
    (phi > 0 ? halfPi : -halfPi) * (by + di * (cy - ay) / 2 + di * di * (cy - 2 * by + ay) / 2)
  ];
}

robinsonRaw.invert = function(x, y) {
  var yy = y / halfPi,
      phi = yy * 90,
      i = min(18, abs(phi / 5)),
      i0 = max(0, floor(i));
  do {
    var ay = robinson_K[i0][1],
        by = robinson_K[i0 + 1][1],
        cy = robinson_K[min(19, i0 + 2)][1],
        u = cy - ay,
        v = cy - 2 * by + ay,
        t = 2 * (abs(yy) - by) / u,
        c = v / u,
        di = t * (1 - c * t * (1 - 2 * c * t));
    if (di >= 0 || i0 === 1) {
      phi = (y >= 0 ? 5 : -5) * (di + i);
      var j = 50, delta;
      do {
        i = min(18, abs(phi) / 5);
        i0 = floor(i);
        di = i - i0;
        ay = robinson_K[i0][1];
        by = robinson_K[i0 + 1][1];
        cy = robinson_K[min(19, i0 + 2)][1];
        phi -= (delta = (y >= 0 ? halfPi : -halfPi) * (by + di * (cy - ay) / 2 + di * di * (cy - 2 * by + ay) / 2) - y) * degrees;
      } while (abs(delta) > epsilon2 && --j > 0);
      break;
    }
  } while (--i0 >= 0);
  var ax = robinson_K[i0][0],
      bx = robinson_K[i0 + 1][0],
      cx = robinson_K[min(19, i0 + 2)][0];
  return [
    x / (bx + di * (cx - ax) / 2 + di * di * (cx - 2 * bx + ax) / 2),
    phi * radians
  ];
};

/* harmony default export */ var robinson = (function() {
  return Object(d3_geo["geoProjection"])(robinsonRaw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/satellite.js



function satelliteVerticalRaw(P) {
  function forward(lambda, phi) {
    var cosPhi = cos(phi),
        k = (P - 1) / (P - cosPhi * cos(lambda));
    return [
      k * cosPhi * sin(lambda),
      k * sin(phi)
    ];
  }

  forward.invert = function(x, y) {
    var rho2 = x * x + y * y,
        rho = sqrt(rho2),
        sinc = (P - sqrt(1 - rho2 * (P + 1) / (P - 1))) / ((P - 1) / rho + rho / (P - 1));
    return [
      atan2(x * sinc, rho * sqrt(1 - sinc * sinc)),
      rho ? asin(y * sinc / rho) : 0
    ];
  };

  return forward;
}

function satelliteRaw(P, omega) {
  var vertical = satelliteVerticalRaw(P);
  if (!omega) return vertical;
  var cosOmega = cos(omega),
      sinOmega = sin(omega);

  function forward(lambda, phi) {
    var coordinates = vertical(lambda, phi),
        y = coordinates[1],
        A = y * sinOmega / (P - 1) + cosOmega;
    return [
      coordinates[0] * cosOmega / A,
      y / A
    ];
  }

  forward.invert = function(x, y) {
    var k = (P - 1) / (P - 1 - y * sinOmega);
    return vertical.invert(k * x, k * y * cosOmega);
  };

  return forward;
}

/* harmony default export */ var satellite = (function() {
  var distance = 2,
      omega = 0,
      m = Object(d3_geo["geoProjectionMutator"])(satelliteRaw),
      p = m(distance, omega);

  // As a multiple of radius.
  p.distance = function(_) {
    if (!arguments.length) return distance;
    return m(distance = +_, omega);
  };

  p.tilt = function(_) {
    if (!arguments.length) return omega * degrees;
    return m(distance, omega = _ * radians);
  };

  return p
      .scale(432.147)
      .clipAngle(acos(1 / distance) * degrees - 1e-6);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/stitch.js
var stitch_epsilon = 1e-4,
    epsilonInverse = 1e4,
    stitch_x0 = -180, x0e = stitch_x0 + stitch_epsilon,
    stitch_x1 = 180, x1e = stitch_x1 - stitch_epsilon,
    stitch_y0 = -90, y0e = stitch_y0 + stitch_epsilon,
    stitch_y1 = 90, y1e = stitch_y1 - stitch_epsilon;

function stitch_quantize(x) {
  return Math.floor(x * epsilonInverse) / epsilonInverse;
}

function normalizePoint(y) {
  return y === stitch_y0 || y === stitch_y1
      ? [0, y] // pole
      : [stitch_x0, stitch_quantize(y)]; // antimeridian
}

function clampPoint(p) {
  if (p[0] <= x0e) p[0] = stitch_x0;
  else if (p[0] >= x1e) p[0] = stitch_x1;
  if (p[1] <= y0e) p[1] = stitch_y0;
  else if (p[1] >= y1e) p[1] = stitch_y1;
}

function clampPoints(points) {
  points.forEach(clampPoint);
}

// For each ring, detect where it crosses the antimeridian or pole.
function extractFragments(polygon, fragments) {
  for (var j = 0, m = polygon.length; j < m; ++j) {
    var ring = polygon[j];
    ring.polygon = polygon;

    // By default, assume that this ring doesn’t need any stitching.
    fragments.push(ring);

    for (var i = 0, n = ring.length; i < n; ++i) {
      var point = ring[i],
          x = point[0],
          y = point[1];

      // If this is an antimeridian or polar point…
      if (x <= x0e || x >= x1e || y <= y0e || y >= y1e) {
        clampPoint(point);

        // Advance through any antimeridian or polar points…
        for (var k = i + 1; k < n; ++k) {
          var pointk = ring[k],
              xk = pointk[0],
              yk = pointk[1];
          if (xk > x0e && xk < x1e && yk > y0e && yk < y1e) break;
        }

        // If this was just a single antimeridian or polar point,
        // we don’t need to cut this ring into a fragment;
        // we can just leave it as-is.
        if (k === i + 1) continue;

        // Otherwise, if this is not the first point in the ring,
        // cut the current fragment so that it ends at the current point.
        // The current point is also normalized for later joining.
        if (i) {
          var fragmentBefore = ring.slice(0, i + 1);
          fragmentBefore.polygon = polygon;
          fragmentBefore[fragmentBefore.length - 1] = normalizePoint(y);
          fragments[fragments.length - 1] = fragmentBefore;
        }

        // If the ring started with an antimeridian fragment,
        // we can ignore that fragment entirely.
        else fragments.pop();

        // If the remainder of the ring is an antimeridian fragment,
        // move on to the next ring.
        if (k >= n) break;

        // Otherwise, add the remaining ring fragment and continue.
        fragments.push(ring = ring.slice(k - 1));
        ring[0] = normalizePoint(ring[0][1]);
        ring.polygon = polygon;
        i = -1;
        n = ring.length;
      }
    }
  }
  polygon.length = 0;
}

// Now stitch the fragments back together into rings.
// TODO remove empty polygons.
function stitchFragments(fragments) {
  var i, n = fragments.length;

  // To connect the fragments start-to-end, create a simple index by end.
  var fragmentByStart = {},
      fragmentByEnd = {},
      fragment,
      start,
      startFragment,
      end,
      endFragment;

  // For each fragment…
  for (i = 0; i < n; ++i) {
    fragment = fragments[i];
    start = fragment[0];
    end = fragment[fragment.length - 1];

    // If this fragment is closed, add it as a standalone ring.
    if (start[0] === end[0] && start[1] === end[1]) {
      fragment.polygon.push(fragment);
      fragments[i] = null;
      continue;
    }

    fragment.index = i;
    fragmentByStart[start] = fragmentByEnd[end] = fragment;
  }

  // For each open fragment…
  for (i = 0; i < n; ++i) {
    fragment = fragments[i];
    if (fragment) {
      start = fragment[0];
      end = fragment[fragment.length - 1];
      startFragment = fragmentByEnd[start];
      endFragment = fragmentByStart[end];

      delete fragmentByStart[start];
      delete fragmentByEnd[end];

      // If this fragment is closed, add it as a standalone ring.
      if (start[0] === end[0] && start[1] === end[1]) {
        fragment.polygon.push(fragment);
        continue;
      }

      if (startFragment) {
        delete fragmentByEnd[start];
        delete fragmentByStart[startFragment[0]];
        startFragment.pop(); // drop the shared coordinate
        fragments[startFragment.index] = null;
        fragment = startFragment.concat(fragment);
        fragment.polygon = startFragment.polygon;

        if (startFragment === endFragment) {
          // Connect both ends to this single fragment to create a ring.
          fragment.polygon.push(fragment);
        } else {
          fragment.index = n++;
          fragments.push(fragmentByStart[fragment[0]] = fragmentByEnd[fragment[fragment.length - 1]] = fragment);
        }
      } else if (endFragment) {
        delete fragmentByStart[end];
        delete fragmentByEnd[endFragment[endFragment.length - 1]];
        fragment.pop(); // drop the shared coordinate
        fragment = fragment.concat(endFragment);
        fragment.polygon = endFragment.polygon;
        fragment.index = n++;
        fragments[endFragment.index] = null;
        fragments.push(fragmentByStart[fragment[0]] = fragmentByEnd[fragment[fragment.length - 1]] = fragment);
      } else {
        fragment.push(fragment[0]); // close ring
        fragment.polygon.push(fragment);
      }
    }
  }
}

function stitchFeature(o) {
  stitchGeometry(o.geometry);
}

function stitchGeometry(o) {
  if (!o) return;
  var fragments, i, n;

  switch (o.type) {
    case "GeometryCollection": {
      o.geometries.forEach(stitchGeometry);
      return;
    }
    case "Point": {
      clampPoint(o.coordinates);
      break;
    }
    case "MultiPoint":
    case "LineString": {
      clampPoints(o.coordinates);
      break;
    }
    case "MultiLineString": {
      o.coordinates.forEach(clampPoints);
      break;
    }
    case "Polygon": {
      extractFragments(o.coordinates, fragments = []);
      break;
    }
    case "MultiPolygon": {
      fragments = [], i = -1, n = o.coordinates.length;
      while (++i < n) extractFragments(o.coordinates[i], fragments);
      break;
    }
    default: return;
  }

  stitchFragments(fragments);
}

/* harmony default export */ var stitch = (function(o) {
  if (o) switch (o.type) {
    case "Feature": stitchFeature(o); break;
    case "FeatureCollection": o.features.forEach(stitchFeature); break;
    default: stitchGeometry(o); break;
  }
  return o;
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/times.js



function timesRaw(lambda, phi) {
  var t = tan(phi / 2),
      s = sin(quarterPi * t);
  return [
    lambda * (0.74482 - 0.34588 * s * s),
    1.70711 * t
  ];
}

timesRaw.invert = function(x, y) {
  var t = y / 1.70711,
      s = sin(quarterPi * t);
  return [
    x / (0.74482 - 0.34588 * s * s),
    2 * atan(t)
  ];
};

/* harmony default export */ var times = (function() {
  return Object(d3_geo["geoProjection"])(timesRaw)
      .scale(146.153);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/twoPoint.js



// Compute the origin as the midpoint of the two reference points.
// Rotate one of the reference points by the origin.
// Apply the spherical law of sines to compute gamma rotation.
/* harmony default export */ var twoPoint = (function(raw, p0, p1) {
  var i = Object(d3_geo["geoInterpolate"])(p0, p1),
      o = i(0.5),
      a = Object(d3_geo["geoRotation"])([-o[0], -o[1]])(p0),
      b = i.distance / 2,
      y = -asin(sin(a[1] * radians) / sin(b)),
      R = [-o[0], -o[1], -(a[0] > 0 ? pi - y : y) * degrees],
      p = Object(d3_geo["geoProjection"])(raw(b)).rotate(R),
      r = Object(d3_geo["geoRotation"])(R),
      center = p.center;

  delete p.rotate;

  p.center = function(_) {
    return arguments.length ? center(r(_)) : r.invert(center());
  };

  return p
      .clipAngle(90);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/twoPointAzimuthal.js




function twoPointAzimuthalRaw(d) {
  var cosd = cos(d);

  function forward(lambda, phi) {
    var coordinates = Object(d3_geo["geoGnomonicRaw"])(lambda, phi);
    coordinates[0] *= cosd;
    return coordinates;
  }

  forward.invert = function(x, y) {
    return d3_geo["geoGnomonicRaw"].invert(x / cosd, y);
  };

  return forward;
}

function twoPointAzimuthalUsa() {
  return twoPointAzimuthal([-158, 21.5], [-77, 39])
      .clipAngle(60)
      .scale(400);
}

function twoPointAzimuthal(p0, p1) {
  return twoPoint(twoPointAzimuthalRaw, p0, p1);
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/twoPointEquidistant.js




// TODO clip to ellipse
function twoPointEquidistantRaw(z0) {
  if (!(z0 *= 2)) return d3_geo["geoAzimuthalEquidistantRaw"];
  var lambdaa = -z0 / 2,
      lambdab = -lambdaa,
      z02 = z0 * z0,
      tanLambda0 = tan(lambdab),
      S = 0.5 / sin(lambdab);

  function forward(lambda, phi) {
    var za = acos(cos(phi) * cos(lambda - lambdaa)),
        zb = acos(cos(phi) * cos(lambda - lambdab)),
        ys = phi < 0 ? -1 : 1;
    za *= za, zb *= zb;
    return [
      (za - zb) / (2 * z0),
      ys * sqrt(4 * z02 * zb - (z02 - za + zb) * (z02 - za + zb)) / (2 * z0)
    ];
  }

  forward.invert = function(x, y) {
    var y2 = y * y,
        cosza = cos(sqrt(y2 + (t = x + lambdaa) * t)),
        coszb = cos(sqrt(y2 + (t = x + lambdab) * t)),
        t,
        d;
    return [
      atan2(d = cosza - coszb, t = (cosza + coszb) * tanLambda0),
      (y < 0 ? -1 : 1) * acos(sqrt(t * t + d * d) * S)
    ];
  };

  return forward;
}

function twoPointEquidistantUsa() {
  return twoPointEquidistant([-158, 21.5], [-77, 39])
      .clipAngle(130)
      .scale(122.571);
}

function twoPointEquidistant(p0, p1) {
  return twoPoint(twoPointEquidistantRaw, p0, p1);
}

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/vanDerGrinten.js



function vanDerGrintenRaw(lambda, phi) {
  if (abs(phi) < math_epsilon) return [lambda, 0];
  var sinTheta = abs(phi / halfPi),
      theta = asin(sinTheta);
  if (abs(lambda) < math_epsilon || abs(abs(phi) - halfPi) < math_epsilon) return [0, sign(phi) * pi * tan(theta / 2)];
  var cosTheta = cos(theta),
      A = abs(pi / lambda - lambda / pi) / 2,
      A2 = A * A,
      G = cosTheta / (sinTheta + cosTheta - 1),
      P = G * (2 / sinTheta - 1),
      P2 = P * P,
      P2_A2 = P2 + A2,
      G_P2 = G - P2,
      Q = A2 + G;
  return [
    sign(lambda) * pi * (A * G_P2 + sqrt(A2 * G_P2 * G_P2 - P2_A2 * (G * G - P2))) / P2_A2,
    sign(phi) * pi * (P * Q - A * sqrt((A2 + 1) * P2_A2 - Q * Q)) / P2_A2
  ];
}

vanDerGrintenRaw.invert = function(x, y) {
  if (abs(y) < math_epsilon) return [x, 0];
  if (abs(x) < math_epsilon) return [0, halfPi * sin(2 * atan(y / pi))];
  var x2 = (x /= pi) * x,
      y2 = (y /= pi) * y,
      x2_y2 = x2 + y2,
      z = x2_y2 * x2_y2,
      c1 = -abs(y) * (1 + x2_y2),
      c2 = c1 - 2 * y2 + x2,
      c3 = -2 * c1 + 1 + 2 * y2 + z,
      d = y2 / c3 + (2 * c2 * c2 * c2 / (c3 * c3 * c3) - 9 * c1 * c2 / (c3 * c3)) / 27,
      a1 = (c1 - c2 * c2 / (3 * c3)) / c3,
      m1 = 2 * sqrt(-a1 / 3),
      theta1 = acos(3 * d / (a1 * m1)) / 3;
  return [
    pi * (x2_y2 - 1 + sqrt(1 + 2 * (x2 - y2) + z)) / (2 * x),
    sign(y) * pi * (-m1 * cos(theta1 + pi / 3) - c2 / (3 * c3))
  ];
};

/* harmony default export */ var vanDerGrinten = (function() {
  return Object(d3_geo["geoProjection"])(vanDerGrintenRaw)
      .scale(79.4183);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/vanDerGrinten2.js



function vanDerGrinten2Raw(lambda, phi) {
  if (abs(phi) < math_epsilon) return [lambda, 0];
  var sinTheta = abs(phi / halfPi),
      theta = asin(sinTheta);
  if (abs(lambda) < math_epsilon || abs(abs(phi) - halfPi) < math_epsilon) return [0, sign(phi) * pi * tan(theta / 2)];
  var cosTheta = cos(theta),
      A = abs(pi / lambda - lambda / pi) / 2,
      A2 = A * A,
      x1 = cosTheta * (sqrt(1 + A2) - A * cosTheta) / (1 + A2 * sinTheta * sinTheta);
  return [
    sign(lambda) * pi * x1,
    sign(phi) * pi * sqrt(1 - x1 * (2 * A + x1))
  ];
}

vanDerGrinten2Raw.invert = function(x, y) {
  if (!x) return [0, halfPi * sin(2 * atan(y / pi))];
  var x1 = abs(x / pi),
      A = (1 - x1 * x1 - (y /= pi) * y) / (2 * x1),
      A2 = A * A,
      B = sqrt(A2 + 1);
  return [
    sign(x) * pi * (B - A),
    sign(y) * halfPi * sin(2 * atan2(sqrt((1 - 2 * A * x1) * (A + B) - x1), sqrt(B + A + x1)))
  ];
};

/* harmony default export */ var vanDerGrinten2 = (function() {
  return Object(d3_geo["geoProjection"])(vanDerGrinten2Raw)
      .scale(79.4183);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/vanDerGrinten3.js



function vanDerGrinten3Raw(lambda, phi) {
  if (abs(phi) < math_epsilon) return [lambda, 0];
  var sinTheta = phi / halfPi,
      theta = asin(sinTheta);
  if (abs(lambda) < math_epsilon || abs(abs(phi) - halfPi) < math_epsilon) return [0, pi * tan(theta / 2)];
  var A = (pi / lambda - lambda / pi) / 2,
      y1 = sinTheta / (1 + cos(theta));
  return [
    pi * (sign(lambda) * sqrt(A * A + 1 - y1 * y1) - A),
    pi * y1
  ];
}

vanDerGrinten3Raw.invert = function(x, y) {
  if (!y) return [x, 0];
  var y1 = y / pi,
      A = (pi * pi * (1 - y1 * y1) - x * x) / (2 * pi * x);
  return [
    x ? pi * (sign(x) * sqrt(A * A + 1) - A) : 0,
    halfPi * sin(2 * atan(y1))
  ];
};

/* harmony default export */ var vanDerGrinten3 = (function() {
  return Object(d3_geo["geoProjection"])(vanDerGrinten3Raw)
        .scale(79.4183);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/vanDerGrinten4.js



function vanDerGrinten4Raw(lambda, phi) {
  if (!phi) return [lambda, 0];
  var phi0 = abs(phi);
  if (!lambda || phi0 === halfPi) return [0, phi];
  var B = phi0 / halfPi,
      B2 = B * B,
      C = (8 * B - B2 * (B2 + 2) - 5) / (2 * B2 * (B - 1)),
      C2 = C * C,
      BC = B * C,
      B_C2 = B2 + C2 + 2 * BC,
      B_3C = B + 3 * C,
      lambda0 = lambda / halfPi,
      lambda1 = lambda0 + 1 / lambda0,
      D = sign(abs(lambda) - halfPi) * sqrt(lambda1 * lambda1 - 4),
      D2 = D * D,
      F = B_C2 * (B2 + C2 * D2 - 1) + (1 - B2) * (B2 * (B_3C * B_3C + 4 * C2) + 12 * BC * C2 + 4 * C2 * C2),
      x1 = (D * (B_C2 + C2 - 1) + 2 * sqrt(F)) / (4 * B_C2 + D2);
  return [
    sign(lambda) * halfPi * x1,
    sign(phi) * halfPi * sqrt(1 + D * abs(x1) - x1 * x1)
  ];
}

vanDerGrinten4Raw.invert = function(x, y) {
  var delta;
  if (!x || !y) return [x, y];
  y /= pi;
  var x1 = sign(x) * x / halfPi,
      D = (x1 * x1 - 1 + 4 * y * y) / abs(x1),
      D2 = D * D,
      B = 2 * y,
      i = 50;
  do {
    var B2 = B * B,
        C = (8 * B - B2 * (B2 + 2) - 5) / (2 * B2 * (B - 1)),
        C_ = (3 * B - B2 * B - 10) / (2 * B2 * B),
        C2 = C * C,
        BC = B * C,
        B_C = B + C,
        B_C2 = B_C * B_C,
        B_3C = B + 3 * C,
        F = B_C2 * (B2 + C2 * D2 - 1) + (1 - B2) * (B2 * (B_3C * B_3C + 4 * C2) + C2 * (12 * BC + 4 * C2)),
        F_ = -2 * B_C * (4 * BC * C2 + (1 - 4 * B2 + 3 * B2 * B2) * (1 + C_) + C2 * (-6 + 14 * B2 - D2 + (-8 + 8 * B2 - 2 * D2) * C_) + BC * (-8 + 12 * B2 + (-10 + 10 * B2 - D2) * C_)),
        sqrtF = sqrt(F),
        f = D * (B_C2 + C2 - 1) + 2 * sqrtF - x1 * (4 * B_C2 + D2),
        f_ = D * (2 * C * C_ + 2 * B_C * (1 + C_)) + F_ / sqrtF - 8 * B_C * (D * (-1 + C2 + B_C2) + 2 * sqrtF) * (1 + C_) / (D2 + 4 * B_C2);
    B -= delta = f / f_;
  } while (delta > math_epsilon && --i > 0);
  return [
    sign(x) * (sqrt(D * D + 4) + D) * pi / 4,
    halfPi * B
  ];
};

/* harmony default export */ var vanDerGrinten4 = (function() {
  return Object(d3_geo["geoProjection"])(vanDerGrinten4Raw)
      .scale(127.16);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/wagner4.js




var wagner4_A = 4 * pi + 3 * sqrt(3),
    wagner4_B = 2 * sqrt(2 * pi * sqrt(3) / wagner4_A);

var wagner4Raw = mollweideBromleyRaw(wagner4_B * sqrt(3) / pi, wagner4_B, wagner4_A / 6);

/* harmony default export */ var wagner4 = (function() {
  return Object(d3_geo["geoProjection"])(wagner4Raw)
      .scale(176.84);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/wagner6.js



function wagner6Raw(lambda, phi) {
  return [lambda * sqrt(1 - 3 * phi * phi / (pi * pi)), phi];
}

wagner6Raw.invert = function(x, y) {
  return [x / sqrt(1 - 3 * y * y / (pi * pi)), y];
};

/* harmony default export */ var wagner6 = (function() {
  return Object(d3_geo["geoProjection"])(wagner6Raw)
      .scale(152.63);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/wagner7.js



function wagner7Raw(lambda, phi) {
  var s = 0.90631 * sin(phi),
      c0 = sqrt(1 - s * s),
      c1 = sqrt(2 / (1 + c0 * cos(lambda /= 3)));
  return [
    2.66723 * c0 * c1 * sin(lambda),
    1.24104 * s * c1
  ];
}

wagner7Raw.invert = function(x, y) {
  var t1 = x / 2.66723,
      t2 = y / 1.24104,
      p = sqrt(t1 * t1 + t2 * t2),
      c = 2 * asin(p / 2);
  return [
    3 * atan2(x * tan(c), 2.66723 * p),
    p && asin(y * sin(c) / (1.24104 * 0.90631 * p))
  ];
};

/* harmony default export */ var wagner7 = (function() {
  return Object(d3_geo["geoProjection"])(wagner7Raw)
      .scale(172.632);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/wiechel.js



function wiechelRaw(lambda, phi) {
  var cosPhi = cos(phi),
      sinPhi = cos(lambda) * cosPhi,
      sin1_Phi = 1 - sinPhi,
      cosLambda = cos(lambda = atan2(sin(lambda) * cosPhi, -sin(phi))),
      sinLambda = sin(lambda);
  cosPhi = sqrt(1 - sinPhi * sinPhi);
  return [
    sinLambda * cosPhi - cosLambda * sin1_Phi,
    -cosLambda * cosPhi - sinLambda * sin1_Phi
  ];
}

wiechelRaw.invert = function(x, y) {
  var w = (x * x + y * y) / -2,
      k = sqrt(-w * (2 + w)),
      b = y * w + x * k,
      a = x * w - y * k,
      D = sqrt(a * a + b * b);
  return [
    atan2(k * b, D * (1 + w)),
    D ? -asin(k * a / D) : 0
  ];
};

/* harmony default export */ var wiechel = (function() {
  return Object(d3_geo["geoProjection"])(wiechelRaw)
      .rotate([0, -90, 45])
      .scale(124.75)
      .clipAngle(180 - 1e-3);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/src/winkel3.js




function winkel3Raw(lambda, phi) {
  var coordinates = aitoffRaw(lambda, phi);
  return [
    (coordinates[0] + lambda / halfPi) / 2,
    (coordinates[1] + phi) / 2
  ];
}

winkel3Raw.invert = function(x, y) {
  var lambda = x, phi = y, i = 25;
  do {
    var cosphi = cos(phi),
        sinphi = sin(phi),
        sin_2phi = sin(2 * phi),
        sin2phi = sinphi * sinphi,
        cos2phi = cosphi * cosphi,
        sinlambda = sin(lambda),
        coslambda_2 = cos(lambda / 2),
        sinlambda_2 = sin(lambda / 2),
        sin2lambda_2 = sinlambda_2 * sinlambda_2,
        C = 1 - cos2phi * coslambda_2 * coslambda_2,
        E = C ? acos(cosphi * coslambda_2) * sqrt(F = 1 / C) : F = 0,
        F,
        fx = 0.5 * (2 * E * cosphi * sinlambda_2 + lambda / halfPi) - x,
        fy = 0.5 * (E * sinphi + phi) - y,
        dxdlambda = 0.5 * F * (cos2phi * sin2lambda_2 + E * cosphi * coslambda_2 * sin2phi) + 0.5 / halfPi,
        dxdphi = F * (sinlambda * sin_2phi / 4 - E * sinphi * sinlambda_2),
        dydlambda = 0.125 * F * (sin_2phi * sinlambda_2 - E * sinphi * cos2phi * sinlambda),
        dydphi = 0.5 * F * (sin2phi * coslambda_2 + E * sin2lambda_2 * cosphi) + 0.5,
        denominator = dxdphi * dydlambda - dydphi * dxdlambda,
        dlambda = (fy * dxdphi - fx * dydphi) / denominator,
        dphi = (fx * dydlambda - fy * dxdlambda) / denominator;
    lambda -= dlambda, phi -= dphi;
  } while ((abs(dlambda) > math_epsilon || abs(dphi) > math_epsilon) && --i > 0);
  return [lambda, phi];
};

/* harmony default export */ var winkel3 = (function() {
  return Object(d3_geo["geoProjection"])(winkel3Raw)
      .scale(158.837);
});

// CONCATENATED MODULE: ./node_modules/d3-geo-projection/index.js
/* concated harmony reexport geoAiry */__webpack_require__.d(__webpack_exports__, "geoAiry", function() { return airy; });
/* concated harmony reexport geoAiryRaw */__webpack_require__.d(__webpack_exports__, "geoAiryRaw", function() { return airyRaw; });
/* concated harmony reexport geoAitoff */__webpack_require__.d(__webpack_exports__, "geoAitoff", function() { return aitoff; });
/* concated harmony reexport geoAitoffRaw */__webpack_require__.d(__webpack_exports__, "geoAitoffRaw", function() { return aitoffRaw; });
/* concated harmony reexport geoArmadillo */__webpack_require__.d(__webpack_exports__, "geoArmadillo", function() { return armadillo; });
/* concated harmony reexport geoArmadilloRaw */__webpack_require__.d(__webpack_exports__, "geoArmadilloRaw", function() { return armadilloRaw; });
/* concated harmony reexport geoAugust */__webpack_require__.d(__webpack_exports__, "geoAugust", function() { return august; });
/* concated harmony reexport geoAugustRaw */__webpack_require__.d(__webpack_exports__, "geoAugustRaw", function() { return augustRaw; });
/* concated harmony reexport geoBaker */__webpack_require__.d(__webpack_exports__, "geoBaker", function() { return baker; });
/* concated harmony reexport geoBakerRaw */__webpack_require__.d(__webpack_exports__, "geoBakerRaw", function() { return bakerRaw; });
/* concated harmony reexport geoBerghaus */__webpack_require__.d(__webpack_exports__, "geoBerghaus", function() { return berghaus; });
/* concated harmony reexport geoBerghausRaw */__webpack_require__.d(__webpack_exports__, "geoBerghausRaw", function() { return berghausRaw; });
/* concated harmony reexport geoBoggs */__webpack_require__.d(__webpack_exports__, "geoBoggs", function() { return boggs; });
/* concated harmony reexport geoBoggsRaw */__webpack_require__.d(__webpack_exports__, "geoBoggsRaw", function() { return boggsRaw; });
/* concated harmony reexport geoBonne */__webpack_require__.d(__webpack_exports__, "geoBonne", function() { return bonne; });
/* concated harmony reexport geoBonneRaw */__webpack_require__.d(__webpack_exports__, "geoBonneRaw", function() { return bonneRaw; });
/* concated harmony reexport geoBottomley */__webpack_require__.d(__webpack_exports__, "geoBottomley", function() { return bottomley; });
/* concated harmony reexport geoBottomleyRaw */__webpack_require__.d(__webpack_exports__, "geoBottomleyRaw", function() { return bottomleyRaw; });
/* concated harmony reexport geoBromley */__webpack_require__.d(__webpack_exports__, "geoBromley", function() { return bromley; });
/* concated harmony reexport geoBromleyRaw */__webpack_require__.d(__webpack_exports__, "geoBromleyRaw", function() { return bromleyRaw; });
/* concated harmony reexport geoChamberlin */__webpack_require__.d(__webpack_exports__, "geoChamberlin", function() { return chamberlin; });
/* concated harmony reexport geoChamberlinRaw */__webpack_require__.d(__webpack_exports__, "geoChamberlinRaw", function() { return chamberlinRaw; });
/* concated harmony reexport geoChamberlinAfrica */__webpack_require__.d(__webpack_exports__, "geoChamberlinAfrica", function() { return chamberlinAfrica; });
/* concated harmony reexport geoCollignon */__webpack_require__.d(__webpack_exports__, "geoCollignon", function() { return collignon; });
/* concated harmony reexport geoCollignonRaw */__webpack_require__.d(__webpack_exports__, "geoCollignonRaw", function() { return collignonRaw; });
/* concated harmony reexport geoCraig */__webpack_require__.d(__webpack_exports__, "geoCraig", function() { return craig; });
/* concated harmony reexport geoCraigRaw */__webpack_require__.d(__webpack_exports__, "geoCraigRaw", function() { return craigRaw; });
/* concated harmony reexport geoCraster */__webpack_require__.d(__webpack_exports__, "geoCraster", function() { return craster; });
/* concated harmony reexport geoCrasterRaw */__webpack_require__.d(__webpack_exports__, "geoCrasterRaw", function() { return crasterRaw; });
/* concated harmony reexport geoCylindricalEqualArea */__webpack_require__.d(__webpack_exports__, "geoCylindricalEqualArea", function() { return cylindricalEqualArea; });
/* concated harmony reexport geoCylindricalEqualAreaRaw */__webpack_require__.d(__webpack_exports__, "geoCylindricalEqualAreaRaw", function() { return cylindricalEqualAreaRaw; });
/* concated harmony reexport geoCylindricalStereographic */__webpack_require__.d(__webpack_exports__, "geoCylindricalStereographic", function() { return cylindricalStereographic; });
/* concated harmony reexport geoCylindricalStereographicRaw */__webpack_require__.d(__webpack_exports__, "geoCylindricalStereographicRaw", function() { return cylindricalStereographicRaw; });
/* concated harmony reexport geoEckert1 */__webpack_require__.d(__webpack_exports__, "geoEckert1", function() { return eckert1; });
/* concated harmony reexport geoEckert1Raw */__webpack_require__.d(__webpack_exports__, "geoEckert1Raw", function() { return eckert1Raw; });
/* concated harmony reexport geoEckert2 */__webpack_require__.d(__webpack_exports__, "geoEckert2", function() { return eckert2; });
/* concated harmony reexport geoEckert2Raw */__webpack_require__.d(__webpack_exports__, "geoEckert2Raw", function() { return eckert2Raw; });
/* concated harmony reexport geoEckert3 */__webpack_require__.d(__webpack_exports__, "geoEckert3", function() { return eckert3; });
/* concated harmony reexport geoEckert3Raw */__webpack_require__.d(__webpack_exports__, "geoEckert3Raw", function() { return eckert3Raw; });
/* concated harmony reexport geoEckert4 */__webpack_require__.d(__webpack_exports__, "geoEckert4", function() { return eckert4; });
/* concated harmony reexport geoEckert4Raw */__webpack_require__.d(__webpack_exports__, "geoEckert4Raw", function() { return eckert4Raw; });
/* concated harmony reexport geoEckert5 */__webpack_require__.d(__webpack_exports__, "geoEckert5", function() { return eckert5; });
/* concated harmony reexport geoEckert5Raw */__webpack_require__.d(__webpack_exports__, "geoEckert5Raw", function() { return eckert5Raw; });
/* concated harmony reexport geoEckert6 */__webpack_require__.d(__webpack_exports__, "geoEckert6", function() { return eckert6; });
/* concated harmony reexport geoEckert6Raw */__webpack_require__.d(__webpack_exports__, "geoEckert6Raw", function() { return eckert6Raw; });
/* concated harmony reexport geoEisenlohr */__webpack_require__.d(__webpack_exports__, "geoEisenlohr", function() { return eisenlohr; });
/* concated harmony reexport geoEisenlohrRaw */__webpack_require__.d(__webpack_exports__, "geoEisenlohrRaw", function() { return eisenlohrRaw; });
/* concated harmony reexport geoFahey */__webpack_require__.d(__webpack_exports__, "geoFahey", function() { return fahey; });
/* concated harmony reexport geoFaheyRaw */__webpack_require__.d(__webpack_exports__, "geoFaheyRaw", function() { return faheyRaw; });
/* concated harmony reexport geoFoucaut */__webpack_require__.d(__webpack_exports__, "geoFoucaut", function() { return foucaut; });
/* concated harmony reexport geoFoucautRaw */__webpack_require__.d(__webpack_exports__, "geoFoucautRaw", function() { return foucautRaw; });
/* concated harmony reexport geoGilbert */__webpack_require__.d(__webpack_exports__, "geoGilbert", function() { return src_gilbert; });
/* concated harmony reexport geoGingery */__webpack_require__.d(__webpack_exports__, "geoGingery", function() { return gingery; });
/* concated harmony reexport geoGingeryRaw */__webpack_require__.d(__webpack_exports__, "geoGingeryRaw", function() { return gingeryRaw; });
/* concated harmony reexport geoGinzburg4 */__webpack_require__.d(__webpack_exports__, "geoGinzburg4", function() { return ginzburg4; });
/* concated harmony reexport geoGinzburg4Raw */__webpack_require__.d(__webpack_exports__, "geoGinzburg4Raw", function() { return ginzburg4Raw; });
/* concated harmony reexport geoGinzburg5 */__webpack_require__.d(__webpack_exports__, "geoGinzburg5", function() { return ginzburg5; });
/* concated harmony reexport geoGinzburg5Raw */__webpack_require__.d(__webpack_exports__, "geoGinzburg5Raw", function() { return ginzburg5Raw; });
/* concated harmony reexport geoGinzburg6 */__webpack_require__.d(__webpack_exports__, "geoGinzburg6", function() { return ginzburg6; });
/* concated harmony reexport geoGinzburg6Raw */__webpack_require__.d(__webpack_exports__, "geoGinzburg6Raw", function() { return ginzburg6Raw; });
/* concated harmony reexport geoGinzburg8 */__webpack_require__.d(__webpack_exports__, "geoGinzburg8", function() { return ginzburg8; });
/* concated harmony reexport geoGinzburg8Raw */__webpack_require__.d(__webpack_exports__, "geoGinzburg8Raw", function() { return ginzburg8Raw; });
/* concated harmony reexport geoGinzburg9 */__webpack_require__.d(__webpack_exports__, "geoGinzburg9", function() { return ginzburg9; });
/* concated harmony reexport geoGinzburg9Raw */__webpack_require__.d(__webpack_exports__, "geoGinzburg9Raw", function() { return ginzburg9Raw; });
/* concated harmony reexport geoGringorten */__webpack_require__.d(__webpack_exports__, "geoGringorten", function() { return gringorten; });
/* concated harmony reexport geoGringortenRaw */__webpack_require__.d(__webpack_exports__, "geoGringortenRaw", function() { return gringortenRaw; });
/* concated harmony reexport geoGuyou */__webpack_require__.d(__webpack_exports__, "geoGuyou", function() { return guyou; });
/* concated harmony reexport geoGuyouRaw */__webpack_require__.d(__webpack_exports__, "geoGuyouRaw", function() { return guyouRaw; });
/* concated harmony reexport geoHammer */__webpack_require__.d(__webpack_exports__, "geoHammer", function() { return hammer; });
/* concated harmony reexport geoHammerRaw */__webpack_require__.d(__webpack_exports__, "geoHammerRaw", function() { return hammerRaw; });
/* concated harmony reexport geoHammerRetroazimuthal */__webpack_require__.d(__webpack_exports__, "geoHammerRetroazimuthal", function() { return hammerRetroazimuthal; });
/* concated harmony reexport geoHammerRetroazimuthalRaw */__webpack_require__.d(__webpack_exports__, "geoHammerRetroazimuthalRaw", function() { return hammerRetroazimuthalRaw; });
/* concated harmony reexport geoHealpix */__webpack_require__.d(__webpack_exports__, "geoHealpix", function() { return healpix; });
/* concated harmony reexport geoHealpixRaw */__webpack_require__.d(__webpack_exports__, "geoHealpixRaw", function() { return healpixRaw; });
/* concated harmony reexport geoHill */__webpack_require__.d(__webpack_exports__, "geoHill", function() { return hill; });
/* concated harmony reexport geoHillRaw */__webpack_require__.d(__webpack_exports__, "geoHillRaw", function() { return hillRaw; });
/* concated harmony reexport geoHomolosine */__webpack_require__.d(__webpack_exports__, "geoHomolosine", function() { return homolosine; });
/* concated harmony reexport geoHomolosineRaw */__webpack_require__.d(__webpack_exports__, "geoHomolosineRaw", function() { return homolosineRaw; });
/* concated harmony reexport geoInterrupt */__webpack_require__.d(__webpack_exports__, "geoInterrupt", function() { return interrupted; });
/* concated harmony reexport geoInterruptedBoggs */__webpack_require__.d(__webpack_exports__, "geoInterruptedBoggs", function() { return interrupted_boggs; });
/* concated harmony reexport geoInterruptedHomolosine */__webpack_require__.d(__webpack_exports__, "geoInterruptedHomolosine", function() { return interrupted_homolosine; });
/* concated harmony reexport geoInterruptedMollweide */__webpack_require__.d(__webpack_exports__, "geoInterruptedMollweide", function() { return interrupted_mollweide; });
/* concated harmony reexport geoInterruptedMollweideHemispheres */__webpack_require__.d(__webpack_exports__, "geoInterruptedMollweideHemispheres", function() { return mollweideHemispheres; });
/* concated harmony reexport geoInterruptedSinuMollweide */__webpack_require__.d(__webpack_exports__, "geoInterruptedSinuMollweide", function() { return interrupted_sinuMollweide; });
/* concated harmony reexport geoInterruptedSinusoidal */__webpack_require__.d(__webpack_exports__, "geoInterruptedSinusoidal", function() { return interrupted_sinusoidal; });
/* concated harmony reexport geoKavrayskiy7 */__webpack_require__.d(__webpack_exports__, "geoKavrayskiy7", function() { return kavrayskiy7; });
/* concated harmony reexport geoKavrayskiy7Raw */__webpack_require__.d(__webpack_exports__, "geoKavrayskiy7Raw", function() { return kavrayskiy7Raw; });
/* concated harmony reexport geoLagrange */__webpack_require__.d(__webpack_exports__, "geoLagrange", function() { return lagrange; });
/* concated harmony reexport geoLagrangeRaw */__webpack_require__.d(__webpack_exports__, "geoLagrangeRaw", function() { return lagrangeRaw; });
/* concated harmony reexport geoLarrivee */__webpack_require__.d(__webpack_exports__, "geoLarrivee", function() { return larrivee; });
/* concated harmony reexport geoLarriveeRaw */__webpack_require__.d(__webpack_exports__, "geoLarriveeRaw", function() { return larriveeRaw; });
/* concated harmony reexport geoLaskowski */__webpack_require__.d(__webpack_exports__, "geoLaskowski", function() { return laskowski; });
/* concated harmony reexport geoLaskowskiRaw */__webpack_require__.d(__webpack_exports__, "geoLaskowskiRaw", function() { return laskowskiRaw; });
/* concated harmony reexport geoLittrow */__webpack_require__.d(__webpack_exports__, "geoLittrow", function() { return littrow; });
/* concated harmony reexport geoLittrowRaw */__webpack_require__.d(__webpack_exports__, "geoLittrowRaw", function() { return littrowRaw; });
/* concated harmony reexport geoLoximuthal */__webpack_require__.d(__webpack_exports__, "geoLoximuthal", function() { return loximuthal; });
/* concated harmony reexport geoLoximuthalRaw */__webpack_require__.d(__webpack_exports__, "geoLoximuthalRaw", function() { return loximuthalRaw; });
/* concated harmony reexport geoMiller */__webpack_require__.d(__webpack_exports__, "geoMiller", function() { return miller; });
/* concated harmony reexport geoMillerRaw */__webpack_require__.d(__webpack_exports__, "geoMillerRaw", function() { return millerRaw; });
/* concated harmony reexport geoModifiedStereographic */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographic", function() { return modifiedStereographic; });
/* concated harmony reexport geoModifiedStereographicRaw */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographicRaw", function() { return modifiedStereographicRaw; });
/* concated harmony reexport geoModifiedStereographicAlaska */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographicAlaska", function() { return modifiedStereographicAlaska; });
/* concated harmony reexport geoModifiedStereographicGs48 */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographicGs48", function() { return modifiedStereographicGs48; });
/* concated harmony reexport geoModifiedStereographicGs50 */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographicGs50", function() { return modifiedStereographicGs50; });
/* concated harmony reexport geoModifiedStereographicMiller */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographicMiller", function() { return modifiedStereographicMiller; });
/* concated harmony reexport geoModifiedStereographicLee */__webpack_require__.d(__webpack_exports__, "geoModifiedStereographicLee", function() { return modifiedStereographicLee; });
/* concated harmony reexport geoMollweide */__webpack_require__.d(__webpack_exports__, "geoMollweide", function() { return mollweide; });
/* concated harmony reexport geoMollweideRaw */__webpack_require__.d(__webpack_exports__, "geoMollweideRaw", function() { return mollweideRaw; });
/* concated harmony reexport geoMtFlatPolarParabolic */__webpack_require__.d(__webpack_exports__, "geoMtFlatPolarParabolic", function() { return mtFlatPolarParabolic; });
/* concated harmony reexport geoMtFlatPolarParabolicRaw */__webpack_require__.d(__webpack_exports__, "geoMtFlatPolarParabolicRaw", function() { return mtFlatPolarParabolicRaw; });
/* concated harmony reexport geoMtFlatPolarQuartic */__webpack_require__.d(__webpack_exports__, "geoMtFlatPolarQuartic", function() { return mtFlatPolarQuartic; });
/* concated harmony reexport geoMtFlatPolarQuarticRaw */__webpack_require__.d(__webpack_exports__, "geoMtFlatPolarQuarticRaw", function() { return mtFlatPolarQuarticRaw; });
/* concated harmony reexport geoMtFlatPolarSinusoidal */__webpack_require__.d(__webpack_exports__, "geoMtFlatPolarSinusoidal", function() { return mtFlatPolarSinusoidal; });
/* concated harmony reexport geoMtFlatPolarSinusoidalRaw */__webpack_require__.d(__webpack_exports__, "geoMtFlatPolarSinusoidalRaw", function() { return mtFlatPolarSinusoidalRaw; });
/* concated harmony reexport geoNaturalEarth */__webpack_require__.d(__webpack_exports__, "geoNaturalEarth", function() { return naturalEarth; });
/* concated harmony reexport geoNaturalEarthRaw */__webpack_require__.d(__webpack_exports__, "geoNaturalEarthRaw", function() { return naturalEarthRaw; });
/* concated harmony reexport geoNellHammer */__webpack_require__.d(__webpack_exports__, "geoNellHammer", function() { return nellHammer; });
/* concated harmony reexport geoNellHammerRaw */__webpack_require__.d(__webpack_exports__, "geoNellHammerRaw", function() { return nellHammerRaw; });
/* concated harmony reexport geoPatterson */__webpack_require__.d(__webpack_exports__, "geoPatterson", function() { return patterson; });
/* concated harmony reexport geoPattersonRaw */__webpack_require__.d(__webpack_exports__, "geoPattersonRaw", function() { return pattersonRaw; });
/* concated harmony reexport geoPolyconic */__webpack_require__.d(__webpack_exports__, "geoPolyconic", function() { return polyconic; });
/* concated harmony reexport geoPolyconicRaw */__webpack_require__.d(__webpack_exports__, "geoPolyconicRaw", function() { return polyconicRaw; });
/* concated harmony reexport geoPolyhedral */__webpack_require__.d(__webpack_exports__, "geoPolyhedral", function() { return polyhedral; });
/* concated harmony reexport geoPolyhedralButterfly */__webpack_require__.d(__webpack_exports__, "geoPolyhedralButterfly", function() { return butterfly; });
/* concated harmony reexport geoPolyhedralCollignon */__webpack_require__.d(__webpack_exports__, "geoPolyhedralCollignon", function() { return polyhedral_collignon; });
/* concated harmony reexport geoPolyhedralWaterman */__webpack_require__.d(__webpack_exports__, "geoPolyhedralWaterman", function() { return waterman; });
/* concated harmony reexport geoProject */__webpack_require__.d(__webpack_exports__, "geoProject", function() { return src_project; });
/* concated harmony reexport geoGringortenQuincuncial */__webpack_require__.d(__webpack_exports__, "geoGringortenQuincuncial", function() { return quincuncial_gringorten; });
/* concated harmony reexport geoPeirceQuincuncial */__webpack_require__.d(__webpack_exports__, "geoPeirceQuincuncial", function() { return peirce; });
/* concated harmony reexport geoPierceQuincuncial */__webpack_require__.d(__webpack_exports__, "geoPierceQuincuncial", function() { return peirce; });
/* concated harmony reexport geoQuantize */__webpack_require__.d(__webpack_exports__, "geoQuantize", function() { return quantize; });
/* concated harmony reexport geoQuincuncial */__webpack_require__.d(__webpack_exports__, "geoQuincuncial", function() { return quincuncial; });
/* concated harmony reexport geoRectangularPolyconic */__webpack_require__.d(__webpack_exports__, "geoRectangularPolyconic", function() { return rectangularPolyconic; });
/* concated harmony reexport geoRectangularPolyconicRaw */__webpack_require__.d(__webpack_exports__, "geoRectangularPolyconicRaw", function() { return rectangularPolyconicRaw; });
/* concated harmony reexport geoRobinson */__webpack_require__.d(__webpack_exports__, "geoRobinson", function() { return robinson; });
/* concated harmony reexport geoRobinsonRaw */__webpack_require__.d(__webpack_exports__, "geoRobinsonRaw", function() { return robinsonRaw; });
/* concated harmony reexport geoSatellite */__webpack_require__.d(__webpack_exports__, "geoSatellite", function() { return satellite; });
/* concated harmony reexport geoSatelliteRaw */__webpack_require__.d(__webpack_exports__, "geoSatelliteRaw", function() { return satelliteRaw; });
/* concated harmony reexport geoSinuMollweide */__webpack_require__.d(__webpack_exports__, "geoSinuMollweide", function() { return sinuMollweide; });
/* concated harmony reexport geoSinuMollweideRaw */__webpack_require__.d(__webpack_exports__, "geoSinuMollweideRaw", function() { return sinuMollweideRaw; });
/* concated harmony reexport geoSinusoidal */__webpack_require__.d(__webpack_exports__, "geoSinusoidal", function() { return sinusoidal; });
/* concated harmony reexport geoSinusoidalRaw */__webpack_require__.d(__webpack_exports__, "geoSinusoidalRaw", function() { return sinusoidalRaw; });
/* concated harmony reexport geoStitch */__webpack_require__.d(__webpack_exports__, "geoStitch", function() { return stitch; });
/* concated harmony reexport geoTimes */__webpack_require__.d(__webpack_exports__, "geoTimes", function() { return times; });
/* concated harmony reexport geoTimesRaw */__webpack_require__.d(__webpack_exports__, "geoTimesRaw", function() { return timesRaw; });
/* concated harmony reexport geoTwoPointAzimuthal */__webpack_require__.d(__webpack_exports__, "geoTwoPointAzimuthal", function() { return twoPointAzimuthal; });
/* concated harmony reexport geoTwoPointAzimuthalRaw */__webpack_require__.d(__webpack_exports__, "geoTwoPointAzimuthalRaw", function() { return twoPointAzimuthalRaw; });
/* concated harmony reexport geoTwoPointAzimuthalUsa */__webpack_require__.d(__webpack_exports__, "geoTwoPointAzimuthalUsa", function() { return twoPointAzimuthalUsa; });
/* concated harmony reexport geoTwoPointEquidistant */__webpack_require__.d(__webpack_exports__, "geoTwoPointEquidistant", function() { return twoPointEquidistant; });
/* concated harmony reexport geoTwoPointEquidistantRaw */__webpack_require__.d(__webpack_exports__, "geoTwoPointEquidistantRaw", function() { return twoPointEquidistantRaw; });
/* concated harmony reexport geoTwoPointEquidistantUsa */__webpack_require__.d(__webpack_exports__, "geoTwoPointEquidistantUsa", function() { return twoPointEquidistantUsa; });
/* concated harmony reexport geoVanDerGrinten */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten", function() { return vanDerGrinten; });
/* concated harmony reexport geoVanDerGrintenRaw */__webpack_require__.d(__webpack_exports__, "geoVanDerGrintenRaw", function() { return vanDerGrintenRaw; });
/* concated harmony reexport geoVanDerGrinten2 */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten2", function() { return vanDerGrinten2; });
/* concated harmony reexport geoVanDerGrinten2Raw */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten2Raw", function() { return vanDerGrinten2Raw; });
/* concated harmony reexport geoVanDerGrinten3 */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten3", function() { return vanDerGrinten3; });
/* concated harmony reexport geoVanDerGrinten3Raw */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten3Raw", function() { return vanDerGrinten3Raw; });
/* concated harmony reexport geoVanDerGrinten4 */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten4", function() { return vanDerGrinten4; });
/* concated harmony reexport geoVanDerGrinten4Raw */__webpack_require__.d(__webpack_exports__, "geoVanDerGrinten4Raw", function() { return vanDerGrinten4Raw; });
/* concated harmony reexport geoWagner4 */__webpack_require__.d(__webpack_exports__, "geoWagner4", function() { return wagner4; });
/* concated harmony reexport geoWagner4Raw */__webpack_require__.d(__webpack_exports__, "geoWagner4Raw", function() { return wagner4Raw; });
/* concated harmony reexport geoWagner6 */__webpack_require__.d(__webpack_exports__, "geoWagner6", function() { return wagner6; });
/* concated harmony reexport geoWagner6Raw */__webpack_require__.d(__webpack_exports__, "geoWagner6Raw", function() { return wagner6Raw; });
/* concated harmony reexport geoWagner7 */__webpack_require__.d(__webpack_exports__, "geoWagner7", function() { return wagner7; });
/* concated harmony reexport geoWagner7Raw */__webpack_require__.d(__webpack_exports__, "geoWagner7Raw", function() { return wagner7Raw; });
/* concated harmony reexport geoWiechel */__webpack_require__.d(__webpack_exports__, "geoWiechel", function() { return wiechel; });
/* concated harmony reexport geoWiechelRaw */__webpack_require__.d(__webpack_exports__, "geoWiechelRaw", function() { return wiechelRaw; });
/* concated harmony reexport geoWinkel3 */__webpack_require__.d(__webpack_exports__, "geoWinkel3", function() { return winkel3; });
/* concated harmony reexport geoWinkel3Raw */__webpack_require__.d(__webpack_exports__, "geoWinkel3Raw", function() { return winkel3Raw; });





































































 // DEPRECATED misspelling






















/***/ }),
/* 22 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
__webpack_require__.r(__webpack_exports__);

// CONCATENATED MODULE: ./node_modules/topojson-client/src/identity.js
/* harmony default export */ var identity = (function(x) {
  return x;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/transform.js


/* harmony default export */ var src_transform = (function(topology) {
  if ((transform = topology.transform) == null) return identity;
  var transform,
      x0,
      y0,
      kx = transform.scale[0],
      ky = transform.scale[1],
      dx = transform.translate[0],
      dy = transform.translate[1];
  return function(point, i) {
    if (!i) x0 = y0 = 0;
    point[0] = (x0 += point[0]) * kx + dx;
    point[1] = (y0 += point[1]) * ky + dy;
    return point;
  };
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/bbox.js


/* harmony default export */ var src_bbox = (function(topology) {
  var bbox = topology.bbox;

  function bboxPoint(p0) {
    p1[0] = p0[0], p1[1] = p0[1], t(p1);
    if (p1[0] < x0) x0 = p1[0];
    if (p1[0] > x1) x1 = p1[0];
    if (p1[1] < y0) y0 = p1[1];
    if (p1[1] > y1) y1 = p1[1];
  }

  function bboxGeometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(bboxGeometry); break;
      case "Point": bboxPoint(o.coordinates); break;
      case "MultiPoint": o.coordinates.forEach(bboxPoint); break;
    }
  }

  if (!bbox) {
    var t = src_transform(topology), p0, p1 = new Array(2), name,
        x0 = Infinity, y0 = x0, x1 = -x0, y1 = -x0;

    topology.arcs.forEach(function(arc) {
      var i = -1, n = arc.length;
      while (++i < n) {
        p0 = arc[i], p1[0] = p0[0], p1[1] = p0[1], t(p1, i);
        if (p1[0] < x0) x0 = p1[0];
        if (p1[0] > x1) x1 = p1[0];
        if (p1[1] < y0) y0 = p1[1];
        if (p1[1] > y1) y1 = p1[1];
      }
    });

    for (name in topology.objects) {
      bboxGeometry(topology.objects[name]);
    }

    bbox = topology.bbox = [x0, y0, x1, y1];
  }

  return bbox;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/reverse.js
/* harmony default export */ var reverse = (function(array, n) {
  var t, j = array.length, i = j - n;
  while (i < --j) t = array[i], array[i++] = array[j], array[j] = t;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/feature.js



/* harmony default export */ var feature = (function(topology, o) {
  return o.type === "GeometryCollection"
      ? {type: "FeatureCollection", features: o.geometries.map(function(o) { return feature_feature(topology, o); })}
      : feature_feature(topology, o);
});

function feature_feature(topology, o) {
  var id = o.id,
      bbox = o.bbox,
      properties = o.properties == null ? {} : o.properties,
      geometry = feature_object(topology, o);
  return id == null && bbox == null ? {type: "Feature", properties: properties, geometry: geometry}
      : bbox == null ? {type: "Feature", id: id, properties: properties, geometry: geometry}
      : {type: "Feature", id: id, bbox: bbox, properties: properties, geometry: geometry};
}

function feature_object(topology, o) {
  var transformPoint = src_transform(topology),
      arcs = topology.arcs;

  function arc(i, points) {
    if (points.length) points.pop();
    for (var a = arcs[i < 0 ? ~i : i], k = 0, n = a.length; k < n; ++k) {
      points.push(transformPoint(a[k].slice(), k));
    }
    if (i < 0) reverse(points, n);
  }

  function point(p) {
    return transformPoint(p.slice());
  }

  function line(arcs) {
    var points = [];
    for (var i = 0, n = arcs.length; i < n; ++i) arc(arcs[i], points);
    if (points.length < 2) points.push(points[0].slice());
    return points;
  }

  function ring(arcs) {
    var points = line(arcs);
    while (points.length < 4) points.push(points[0].slice());
    return points;
  }

  function polygon(arcs) {
    return arcs.map(ring);
  }

  function geometry(o) {
    var type = o.type, coordinates;
    switch (type) {
      case "GeometryCollection": return {type: type, geometries: o.geometries.map(geometry)};
      case "Point": coordinates = point(o.coordinates); break;
      case "MultiPoint": coordinates = o.coordinates.map(point); break;
      case "LineString": coordinates = line(o.arcs); break;
      case "MultiLineString": coordinates = o.arcs.map(line); break;
      case "Polygon": coordinates = polygon(o.arcs); break;
      case "MultiPolygon": coordinates = o.arcs.map(polygon); break;
      default: return null;
    }
    return {type: type, coordinates: coordinates};
  }

  return geometry(o);
}

// CONCATENATED MODULE: ./node_modules/topojson-client/src/stitch.js
/* harmony default export */ var stitch = (function(topology, arcs) {
  var stitchedArcs = {},
      fragmentByStart = {},
      fragmentByEnd = {},
      fragments = [],
      emptyIndex = -1;

  // Stitch empty arcs first, since they may be subsumed by other arcs.
  arcs.forEach(function(i, j) {
    var arc = topology.arcs[i < 0 ? ~i : i], t;
    if (arc.length < 3 && !arc[1][0] && !arc[1][1]) {
      t = arcs[++emptyIndex], arcs[emptyIndex] = i, arcs[j] = t;
    }
  });

  arcs.forEach(function(i) {
    var e = ends(i),
        start = e[0],
        end = e[1],
        f, g;

    if (f = fragmentByEnd[start]) {
      delete fragmentByEnd[f.end];
      f.push(i);
      f.end = end;
      if (g = fragmentByStart[end]) {
        delete fragmentByStart[g.start];
        var fg = g === f ? f : f.concat(g);
        fragmentByStart[fg.start = f.start] = fragmentByEnd[fg.end = g.end] = fg;
      } else {
        fragmentByStart[f.start] = fragmentByEnd[f.end] = f;
      }
    } else if (f = fragmentByStart[end]) {
      delete fragmentByStart[f.start];
      f.unshift(i);
      f.start = start;
      if (g = fragmentByEnd[start]) {
        delete fragmentByEnd[g.end];
        var gf = g === f ? f : g.concat(f);
        fragmentByStart[gf.start = g.start] = fragmentByEnd[gf.end = f.end] = gf;
      } else {
        fragmentByStart[f.start] = fragmentByEnd[f.end] = f;
      }
    } else {
      f = [i];
      fragmentByStart[f.start = start] = fragmentByEnd[f.end = end] = f;
    }
  });

  function ends(i) {
    var arc = topology.arcs[i < 0 ? ~i : i], p0 = arc[0], p1;
    if (topology.transform) p1 = [0, 0], arc.forEach(function(dp) { p1[0] += dp[0], p1[1] += dp[1]; });
    else p1 = arc[arc.length - 1];
    return i < 0 ? [p1, p0] : [p0, p1];
  }

  function flush(fragmentByEnd, fragmentByStart) {
    for (var k in fragmentByEnd) {
      var f = fragmentByEnd[k];
      delete fragmentByStart[f.start];
      delete f.start;
      delete f.end;
      f.forEach(function(i) { stitchedArcs[i < 0 ? ~i : i] = 1; });
      fragments.push(f);
    }
  }

  flush(fragmentByEnd, fragmentByStart);
  flush(fragmentByStart, fragmentByEnd);
  arcs.forEach(function(i) { if (!stitchedArcs[i < 0 ? ~i : i]) fragments.push([i]); });

  return fragments;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/mesh.js



/* harmony default export */ var mesh = (function(topology) {
  return feature_object(topology, meshArcs.apply(this, arguments));
});

function meshArcs(topology, object, filter) {
  var arcs, i, n;
  if (arguments.length > 1) arcs = extractArcs(topology, object, filter);
  else for (i = 0, arcs = new Array(n = topology.arcs.length); i < n; ++i) arcs[i] = i;
  return {type: "MultiLineString", arcs: stitch(topology, arcs)};
}

function extractArcs(topology, object, filter) {
  var arcs = [],
      geomsByArc = [],
      geom;

  function extract0(i) {
    var j = i < 0 ? ~i : i;
    (geomsByArc[j] || (geomsByArc[j] = [])).push({i: i, g: geom});
  }

  function extract1(arcs) {
    arcs.forEach(extract0);
  }

  function extract2(arcs) {
    arcs.forEach(extract1);
  }

  function extract3(arcs) {
    arcs.forEach(extract2);
  }

  function geometry(o) {
    switch (geom = o, o.type) {
      case "GeometryCollection": o.geometries.forEach(geometry); break;
      case "LineString": extract1(o.arcs); break;
      case "MultiLineString": case "Polygon": extract2(o.arcs); break;
      case "MultiPolygon": extract3(o.arcs); break;
    }
  }

  geometry(object);

  geomsByArc.forEach(filter == null
      ? function(geoms) { arcs.push(geoms[0].i); }
      : function(geoms) { if (filter(geoms[0].g, geoms[geoms.length - 1].g)) arcs.push(geoms[0].i); });

  return arcs;
}

// CONCATENATED MODULE: ./node_modules/topojson-client/src/merge.js



function planarRingArea(ring) {
  var i = -1, n = ring.length, a, b = ring[n - 1], area = 0;
  while (++i < n) a = b, b = ring[i], area += a[0] * b[1] - a[1] * b[0];
  return Math.abs(area); // Note: doubled area!
}

/* harmony default export */ var merge = (function(topology) {
  return feature_object(topology, mergeArcs.apply(this, arguments));
});

function mergeArcs(topology, objects) {
  var polygonsByArc = {},
      polygons = [],
      groups = [];

  objects.forEach(geometry);

  function geometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(geometry); break;
      case "Polygon": extract(o.arcs); break;
      case "MultiPolygon": o.arcs.forEach(extract); break;
    }
  }

  function extract(polygon) {
    polygon.forEach(function(ring) {
      ring.forEach(function(arc) {
        (polygonsByArc[arc = arc < 0 ? ~arc : arc] || (polygonsByArc[arc] = [])).push(polygon);
      });
    });
    polygons.push(polygon);
  }

  function area(ring) {
    return planarRingArea(feature_object(topology, {type: "Polygon", arcs: [ring]}).coordinates[0]);
  }

  polygons.forEach(function(polygon) {
    if (!polygon._) {
      var group = [],
          neighbors = [polygon];
      polygon._ = 1;
      groups.push(group);
      while (polygon = neighbors.pop()) {
        group.push(polygon);
        polygon.forEach(function(ring) {
          ring.forEach(function(arc) {
            polygonsByArc[arc < 0 ? ~arc : arc].forEach(function(polygon) {
              if (!polygon._) {
                polygon._ = 1;
                neighbors.push(polygon);
              }
            });
          });
        });
      }
    }
  });

  polygons.forEach(function(polygon) {
    delete polygon._;
  });

  return {
    type: "MultiPolygon",
    arcs: groups.map(function(polygons) {
      var arcs = [], n;

      // Extract the exterior (unique) arcs.
      polygons.forEach(function(polygon) {
        polygon.forEach(function(ring) {
          ring.forEach(function(arc) {
            if (polygonsByArc[arc < 0 ? ~arc : arc].length < 2) {
              arcs.push(arc);
            }
          });
        });
      });

      // Stitch the arcs into one or more rings.
      arcs = stitch(topology, arcs);

      // If more than one ring is returned,
      // at most one of these rings can be the exterior;
      // choose the one with the greatest absolute area.
      if ((n = arcs.length) > 1) {
        for (var i = 1, k = area(arcs[0]), ki, t; i < n; ++i) {
          if ((ki = area(arcs[i])) > k) {
            t = arcs[0], arcs[0] = arcs[i], arcs[i] = t, k = ki;
          }
        }
      }

      return arcs;
    })
  };
}

// CONCATENATED MODULE: ./node_modules/topojson-client/src/bisect.js
/* harmony default export */ var bisect = (function(a, x) {
  var lo = 0, hi = a.length;
  while (lo < hi) {
    var mid = lo + hi >>> 1;
    if (a[mid] < x) lo = mid + 1;
    else hi = mid;
  }
  return lo;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/neighbors.js


/* harmony default export */ var src_neighbors = (function(objects) {
  var indexesByArc = {}, // arc index -> array of object indexes
      neighbors = objects.map(function() { return []; });

  function line(arcs, i) {
    arcs.forEach(function(a) {
      if (a < 0) a = ~a;
      var o = indexesByArc[a];
      if (o) o.push(i);
      else indexesByArc[a] = [i];
    });
  }

  function polygon(arcs, i) {
    arcs.forEach(function(arc) { line(arc, i); });
  }

  function geometry(o, i) {
    if (o.type === "GeometryCollection") o.geometries.forEach(function(o) { geometry(o, i); });
    else if (o.type in geometryType) geometryType[o.type](o.arcs, i);
  }

  var geometryType = {
    LineString: line,
    MultiLineString: polygon,
    Polygon: polygon,
    MultiPolygon: function(arcs, i) { arcs.forEach(function(arc) { polygon(arc, i); }); }
  };

  objects.forEach(geometry);

  for (var i in indexesByArc) {
    for (var indexes = indexesByArc[i], m = indexes.length, j = 0; j < m; ++j) {
      for (var k = j + 1; k < m; ++k) {
        var ij = indexes[j], ik = indexes[k], n;
        if ((n = neighbors[ij])[i = bisect(n, ik)] !== ik) n.splice(i, 0, ik);
        if ((n = neighbors[ik])[i = bisect(n, ij)] !== ij) n.splice(i, 0, ij);
      }
    }
  }

  return neighbors;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/quantize.js


/* harmony default export */ var quantize = (function(topology, n) {
  if (!((n = Math.floor(n)) >= 2)) throw new Error("n must be ≥2");
  if (topology.transform) throw new Error("already quantized");
  var bb = src_bbox(topology), name,
      dx = bb[0], kx = (bb[2] - dx) / (n - 1) || 1,
      dy = bb[1], ky = (bb[3] - dy) / (n - 1) || 1;

  function quantizePoint(p) {
    p[0] = Math.round((p[0] - dx) / kx);
    p[1] = Math.round((p[1] - dy) / ky);
  }

  function quantizeGeometry(o) {
    switch (o.type) {
      case "GeometryCollection": o.geometries.forEach(quantizeGeometry); break;
      case "Point": quantizePoint(o.coordinates); break;
      case "MultiPoint": o.coordinates.forEach(quantizePoint); break;
    }
  }

  topology.arcs.forEach(function(arc) {
    var i = 1,
        j = 1,
        n = arc.length,
        pi = arc[0],
        x0 = pi[0] = Math.round((pi[0] - dx) / kx),
        y0 = pi[1] = Math.round((pi[1] - dy) / ky),
        pj,
        x1,
        y1;

    for (; i < n; ++i) {
      pi = arc[i];
      x1 = Math.round((pi[0] - dx) / kx);
      y1 = Math.round((pi[1] - dy) / ky);
      if (x1 !== x0 || y1 !== y0) {
        pj = arc[j++];
        pj[0] = x1 - x0, x0 = x1;
        pj[1] = y1 - y0, y0 = y1;
      }
    }

    if (j < 2) {
      pj = arc[j++];
      pj[0] = 0;
      pj[1] = 0;
    }

    arc.length = j;
  });

  for (name in topology.objects) {
    quantizeGeometry(topology.objects[name]);
  }

  topology.transform = {
    scale: [kx, ky],
    translate: [dx, dy]
  };

  return topology;
});

// CONCATENATED MODULE: ./node_modules/topojson-client/src/untransform.js


/* harmony default export */ var untransform = (function(topology) {
  if ((transform = topology.transform) == null) return identity;
  var transform,
      x0,
      y0,
      kx = transform.scale[0],
      ky = transform.scale[1],
      dx = transform.translate[0],
      dy = transform.translate[1];
  return function(point, i) {
    if (!i) x0 = y0 = 0;
    var x1 = Math.round((point[0] - dx) / kx),
        y1 = Math.round((point[1] - dy) / ky);
    point[0] = x1 - x0, x0 = x1;
    point[1] = y1 - y0, y0 = y1;
    return point;
  };
});

// CONCATENATED MODULE: ./node_modules/topojson-client/index.js
/* concated harmony reexport bbox */__webpack_require__.d(__webpack_exports__, "bbox", function() { return src_bbox; });
/* concated harmony reexport feature */__webpack_require__.d(__webpack_exports__, "feature", function() { return feature; });
/* concated harmony reexport mesh */__webpack_require__.d(__webpack_exports__, "mesh", function() { return mesh; });
/* concated harmony reexport meshArcs */__webpack_require__.d(__webpack_exports__, "meshArcs", function() { return meshArcs; });
/* concated harmony reexport merge */__webpack_require__.d(__webpack_exports__, "merge", function() { return merge; });
/* concated harmony reexport mergeArcs */__webpack_require__.d(__webpack_exports__, "mergeArcs", function() { return mergeArcs; });
/* concated harmony reexport neighbors */__webpack_require__.d(__webpack_exports__, "neighbors", function() { return src_neighbors; });
/* concated harmony reexport quantize */__webpack_require__.d(__webpack_exports__, "quantize", function() { return quantize; });
/* concated harmony reexport transform */__webpack_require__.d(__webpack_exports__, "transform", function() { return src_transform; });
/* concated harmony reexport untransform */__webpack_require__.d(__webpack_exports__, "untransform", function() { return untransform; });










/***/ }),
/* 23 */
/***/ (function(module, __webpack_exports__, __webpack_require__) {

"use strict";
__webpack_require__.r(__webpack_exports__);

// EXTERNAL MODULE: external {"root":"React","commonjs":"react","commonjs2":"react","amd":"react"}
var external_root_React_commonjs_react_commonjs2_react_amd_react_ = __webpack_require__(1);
var external_root_React_commonjs_react_commonjs2_react_amd_react_default = /*#__PURE__*/__webpack_require__.n(external_root_React_commonjs_react_commonjs2_react_amd_react_);

// EXTERNAL MODULE: external {"root":"PropTypes","commonjs":"prop-types","commonjs2":"prop-types","amd":"prop-types"}
var external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_ = __webpack_require__(3);
var external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default = /*#__PURE__*/__webpack_require__.n(external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_);

// EXTERNAL MODULE: ./node_modules/react-simple-maps/lib/index.js
var lib = __webpack_require__(4);

// CONCATENATED MODULE: ./src/map_data.js
/* harmony default export */ var map_data = ([{
  location: "Atlanta",
  lat: 33.748995,
  lng: -84.387982,
  country: "USA",
  chapter: true
}, {
  location: "Austin",
  lat: 30.267153,
  lng: -97.743061,
  country: "USA",
  chapter: true
}, {
  location: "Australia",
  lat: -25.274398,
  lng: 133.775136,
  country: "Australia",
  chapter: true
}, {
  location: "Baltimore",
  lat: 39.290385,
  lng: -76.612189,
  country: "USA",
  chapter: true
}, {
  location: "Bay Area",
  lat: 37.827178,
  lng: -122.291308,
  country: "USA",
  chapter: true
}, {
  location: "Boston",
  lat: 42.360083,
  lng: -71.05888,
  country: "USA",
  chapter: true
}, {
  location: "Brasil",
  lat: -14.235004,
  lng: -51.92528,
  country: "Brazil",
  chapter: true
}, {
  location: "Chicago",
  lat: 41.878114,
  lng: -87.629798,
  country: "USA",
  chapter: true
}, {
  location: "DC",
  lat: 38.907192,
  lng: -77.036871,
  country: "USA",
  chapter: true
}, {
  location: "Dallas",
  lat: 32.776664,
  lng: -96.796988,
  country: "USA",
  chapter: true
}, {
  location: "Detroit",
  lat: 42.331427,
  lng: -83.045754,
  country: "USA",
  chapter: true
}, {
  location: "Deutschland",
  lat: 51.165691,
  lng: 10.451526,
  country: "Germany",
  chapter: true
}, {
  location: "GWU UK",
  lat: 55.378051,
  lng: -3.435973,
  country: "United Kingdom",
  chapter: true
}, {
  location: "Los Angeles",
  lat: 34.052234,
  lng: -118.243685,
  country: "USA",
  chapter: true
}, {
  location: "Montréal",
  lat: 45.501689,
  lng: -73.567256,
  country: "Canada",
  chapter: true
}, {
  location: "New York City",
  lat: 40.712775,
  lng: -74.005973,
  country: "USA",
  chapter: true
}, {
  location: "Orange County",
  lat: 33.717471,
  lng: -117.831143,
  country: "USA",
  chapter: true
}, {
  location: "Orlando",
  lat: 28.538336,
  lng: -81.379236,
  country: "USA",
  chapter: true
}, {
  location: "Ottawa",
  lat: 45.42153,
  lng: -75.697193,
  country: "Canada",
  chapter: true
}, {
  location: "STJV",
  lat: 46.227638,
  lng: 2.213749,
  country: "France",
  chapter: true
}, {
  location: "Seattle",
  lat: 47.60621,
  lng: -122.332071,
  country: "USA",
  chapter: true
}, {
  location: "Spain",
  lat: 40.463667,
  lng: -3.74922,
  country: "Spain",
  chapter: true
}, {
  location: "Sweden",
  lat: 60.128161,
  lng: 18.643501,
  country: "Sweden",
  chapter: true
}, {
  location: "Toronto",
  lat: 43.653226,
  lng: -79.383184,
  country: "Canada",
  chapter: true
}, {
  location: "Triangle",
  lat: 35.899168,
  lng: -78.86364,
  country: "USA",
  chapter: true
}, {
  location: "Vancouver",
  lat: 49.282729,
  lng: -123.120738,
  country: "Canada",
  chapter: true
}, {
  location: "Italy",
  lat: 41.87194,
  lng: 12.56738,
  country: "Italy",
  chapter: false
}, {
  location: "Nashville",
  lat: 36.162664,
  lng: -86.781602,
  country: "USA",
  chapter: false
}, {
  location: "Russia",
  lat: 55.755826,
  lng: 37.6173,
  country: "Russia",
  chapter: true
}, {
  location: "Taipei",
  lat: 25.032969,
  lng: 121.565418,
  country: "Taiwan",
  chapter: false
}, {
  location: "Mexico",
  lat: 23.634501,
  lng: -102.552784,
  country: "Mexico",
  chapter: false
}, {
  location: "Prague",
  lat: 49.817492,
  lng: 15.472962,
  country: "Czech Republic",
  chapter: false
}, {
  location: "Indiana",
  lat: 40.267194,
  lng: -86.134902,
  country: "USA",
  chapter: false
}, {
  location: "Utah",
  lat: 39.32098,
  lng: -111.093731,
  country: "USA",
  chapter: false
}, {
  location: "Bucaramanga",
  lat: 7.119349,
  lng: -73.122742,
  country: "Colombia",
  chapter: false
}, {
  location: "Lewiston, Maine",
  lat: 44.100351,
  lng: -70.214776,
  country: "USA",
  chapter: false
}, {
  location: "Madison, Wisconsin",
  lat: 43.073052,
  lng: -89.40123,
  country: "USA",
  chapter: false
}, {
  location: "Poland",
  lat: 51.919438,
  lng: 19.145136,
  country: "Poland",
  chapter: true
}, {
  location: "Ukraine",
  lat: 48.379433,
  lng: 31.16558,
  country: "Ukraine",
  chapter: true
}, {
  location: "Las Vegas",
  lat: 36.169941,
  lng: -115.13983,
  country: "USA",
  chapter: false
}, {
  location: "Ireland",
  lat: 53.349805,
  lng: -6.26031,
  country: "Ireland",
  chapter: true
}, {
  location: "Alberta",
  lat: 53.933271,
  lng: -116.576503,
  country: "Canada",
  chapter: false
}, {
  location: "Thailand",
  lat: 15.870032,
  lng: 100.992541,
  country: "Thailand",
  chapter: false
}, {
  location: "India",
  lat: 20.593684,
  lng: 78.96288,
  country: "India",
  chapter: false
}, {
  location: "Japan",
  lat: 36.204824,
  lng: 138.252924,
  country: "Japan",
  chapter: false
}, {
  location: "Korea",
  lat: 37.663998,
  lng: 127.978458,
  country: "Korea",
  chapter: false
}, {
  location: "Malaysia",
  lat: 4.210484,
  lng: 101.975766,
  country: "Malaysia",
  chapter: true
}, {
  location: "Vietnam",
  lat: 14.058324,
  lng: 108.277199,
  country: "Vietnam",
  chapter: false
}, {
  location: "The Netherlands",
  lat: 52.132633,
  lng: 5.291266,
  country: "The Netherlands",
  chapter: false
}, {
  location: "Belgium",
  lat: 50.503887,
  lng: 4.469936,
  country: "Belgium",
  chapter: false
}, {
  location: "Finland",
  lat: 61.92411,
  lng: 25.748151,
  country: "Finland",
  chapter: false
}, {
  location: "Argentina",
  lat: -38.416097,
  lng: -63.616672,
  country: "Argentina",
  chapter: true
}, {
  location: "New Zealand",
  lat: -40.900557,
  lng: 174.885971,
  country: "New Zealand",
  chapter: false
}]);
// CONCATENATED MODULE: ./src/GWUMarker.js
function _extends() { _extends = Object.assign || function (target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i]; for (var key in source) { if (Object.prototype.hasOwnProperty.call(source, key)) { target[key] = source[key]; } } } return target; }; return _extends.apply(this, arguments); }

function _objectWithoutPropertiesLoose(source, excluded) { if (source == null) return {}; var target = {}; var sourceKeys = Object.keys(source); var key, i; for (i = 0; i < sourceKeys.length; i++) { key = sourceKeys[i]; if (excluded.indexOf(key) >= 0) continue; target[key] = source[key]; } return target; }

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
        rest = _objectWithoutPropertiesLoose(_this$props, ["marker", "scale"]);

    return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(lib["Marker"], _extends({
      key: marker.name,
      marker: marker
    }, rest), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("g", {
      transform: "scale(" + scale + "), translate(-110, -600)"
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("g", null, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("path", {
      d: "\n                M104.382,145.717c-89.734,0-162.725,73.003-162.725,162.732c0,36.571,25.572,94.673,76.01,172.686\n                c36.714,56.783,73.96,104.287,74.33,104.75c2.985,3.806,7.553,6.017,12.379,6.017c4.829,0,9.396-2.211,12.381-6.017\n                c0.37-0.463,37.622-47.971,74.33-104.75c50.438-78.013,76.012-136.121,76.012-172.686\n                C267.111,218.72,194.117,145.717,104.382,145.717z M42.945,330.743c-3.481,2.971-7.626,5.503-12.439,7.6\n                c-4.804,2.094-10.303,3.146-16.504,3.146c-5.293,0-10.114-0.871-14.47-2.616c-4.357-1.741-8.101-4.16-11.226-7.259\n                c-3.136-3.098-5.571-6.775-7.311-11.037c-1.742-4.259-2.616-8.903-2.616-13.936v-0.197c0-4.838,0.89-9.372,2.662-13.602\n                c1.772-4.222,4.248-7.937,7.405-11.131c3.16-3.193,6.909-5.712,11.234-7.55c4.317-1.839,9.061-2.758,14.228-2.758\n                c3.033,0,5.791,0.208,8.28,0.631c2.481,0.417,4.787,1.016,6.916,1.788c2.136,0.776,4.128,1.745,6.002,2.907\n                c1.876,1.158,3.682,2.481,5.421,3.966l-9.391,11.327c-1.291-1.098-2.583-2.067-3.868-2.906c-1.293-0.834-2.63-1.548-4.018-2.126\n                c-1.395-0.583-2.891-1.033-4.507-1.356c-1.614-0.322-3.387-0.484-5.316-0.484c-2.717,0-5.247,0.565-7.603,1.697\n                c-2.354,1.13-4.419,2.651-6.2,4.562c-1.771,1.906-3.158,4.139-4.158,6.692c-1.001,2.555-1.498,5.289-1.498,8.2v0.193\n                c0,3.104,0.497,5.968,1.498,8.588c1,2.617,2.418,4.898,4.263,6.838c1.835,1.941,4,3.44,6.481,4.513\n                c2.481,1.068,5.247,1.598,8.28,1.598c5.545,0,10.226-1.366,14.03-4.096v-9.749H13.522v-12.872h29.423V330.743z M131.72,340.813\n                h-12.977l-15.29-44.337l-15.299,44.337H75.181l-23.13-68.245H68.02l14.038,45.883l15.195-46.077h12.779l15.195,46.077\n                l14.039-45.883h15.58L131.72,340.813z M221.794,324.256c-1.395,3.81-3.381,6.988-5.956,9.536c-2.585,2.553-5.712,4.455-9.39,5.711\n                c-3.68,1.26-7.777,1.887-12.298,1.887c-9.099,0-16.262-2.517-21.49-7.551c-5.224-5.033-7.839-12.55-7.839-22.554v-38.718h14.905\n                V310.9c0,5.549,1.291,9.729,3.876,12.532c2.575,2.809,6.161,4.217,10.744,4.217c4.579,0,8.161-1.359,10.746-4.07\n                c2.575-2.71,3.868-6.775,3.868-12.194v-38.817h14.912v38.234C223.874,315.965,223.179,320.451,221.794,324.256z\n              "
    }))));
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
  scale: 0.1
};
/* harmony default export */ var src_GWUMarker = (GWUMarker_GWUMarker);
// CONCATENATED MODULE: ./src/ChapterMapComponent.js
function _objectSpread(target) { for (var i = 1; i < arguments.length; i++) { var source = arguments[i] != null ? arguments[i] : {}; var ownKeys = Object.keys(source); if (typeof Object.getOwnPropertySymbols === 'function') { ownKeys = ownKeys.concat(Object.getOwnPropertySymbols(source).filter(function (sym) { return Object.getOwnPropertyDescriptor(source, sym).enumerable; })); } ownKeys.forEach(function (key) { _defineProperty(target, key, source[key]); }); } return target; }

function _defineProperty(obj, key, value) { if (key in obj) { Object.defineProperty(obj, key, { value: value, enumerable: true, configurable: true, writable: true }); } else { obj[key] = value; } return obj; }

function ChapterMapComponent_inheritsLoose(subClass, superClass) { subClass.prototype = Object.create(superClass.prototype); subClass.prototype.constructor = subClass; subClass.__proto__ = superClass; }






var smallWorldDataPromise;
var worldDataPromise;

function getWorldData() {
  return worldDataPromise = worldDataPromise || fetch(location.pathname + "/dist/world-50m.json").then(function (res) {
    return res.json();
  });
}

function getSmallWorldData() {
  return smallWorldDataPromise = smallWorldDataPromise || fetch(location.pathname + "/dist/world-110m.json").then(function (res) {
    return res.json();
  });
}

var markers = map_data.filter(function (point) {
  return point.chapter;
}).map(function (point) {
  return {
    name: point.location,
    coordinates: [point.lng, point.lat]
  };
});
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
      worldData: null
    };
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
  };

  _proto.render = function render() {
    var _this$props = this.props,
        centerLat = _this$props.centerLat,
        centerLng = _this$props.centerLng,
        width = _this$props.width,
        height = _this$props.height,
        scale = _this$props.scale,
        isGeographyIncluded = _this$props.isGeographyIncluded,
        markerScale = _this$props.markerScale,
        style = _this$props.style,
        className = _this$props.className;
    var worldData = this.state.worldData;
    return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("div", {
      className: className,
      style: _objectSpread({}, style || {}, {
        width: width,
        height: height
      })
    }, !worldData && external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement("div", null, "Loading..."), worldData && external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(lib["ComposableMap"], {
      projectionConfig: {
        scale: scale
      },
      width: width,
      height: height,
      style: {
        width: "100%",
        height: "auto"
      }
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(lib["ZoomableGroup"], {
      center: [centerLng, centerLat],
      disablePanning: true
    }, external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(lib["Geographies"], {
      geography: this.state.worldData
    }, function (geographies, projection) {
      return geographies.filter(isGeographyIncluded).map(function (geography, i) {
        var isHighlighted = map_data.some(function (point) {
          return geographyMatchesCountryString(geography, point.country);
        });
        return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(lib["Geography"], {
          key: i,
          geography: geography,
          projection: projection,
          style: {
            default: {
              fill: isHighlighted ? "#ccc" : "#efefef",
              stroke: "#333",
              strokeWidth: 0.3,
              outline: "none"
            },
            hover: {
              fill: "#607D8B",
              stroke: "#607D8B",
              strokeWidth: 0.75,
              outline: "none"
            },
            pressed: {
              fill: "#FF5722",
              stroke: "#607D8B",
              strokeWidth: 0.75,
              outline: "none"
            }
          }
        });
      });
    }), external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(lib["Markers"], null, markers.map(function (marker) {
      return external_root_React_commonjs_react_commonjs2_react_amd_react_default.a.createElement(src_GWUMarker, {
        key: marker.name,
        marker: marker,
        scale: markerScale
      });
    })))));
  };

  return ChapterMapComponent;
}(external_root_React_commonjs_react_commonjs2_react_amd_react_["PureComponent"]);

ChapterMapComponent_ChapterMapComponent.propTypes = {
  centerLat: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  centerLng: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  width: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  height: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  scale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired,
  isGeographyIncluded: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.func.isRequired,
  markerScale: external_root_PropTypes_commonjs_prop_types_commonjs2_prop_types_amd_prop_types_default.a.number.isRequired
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
  markerScale: 0.1
};
/* harmony default export */ var src_ChapterMapComponent = (ChapterMapComponent_ChapterMapComponent);
// CONCATENATED MODULE: ./src/index.js

/* harmony default export */ var src = __webpack_exports__["default"] = (src_ChapterMapComponent);

/***/ })
/******/ ])["default"];
});
//# sourceMappingURL=chapter-map-component.js.map