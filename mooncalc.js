/*
 Copyright (c) 2011-2013, Vladimir Agafonkin
 MoonCalc is a JavaScript library for calculating sun position, sunlight phases, and moon position.
 https://github.com/mourner/MoonCalc
 */

(function (global) { /*jshint smarttabs: true */

	"use strict";


	// export either as a CommonJS module or a global variable

	var MoonCalc;

	if (typeof exports !== 'undefined') {
		MoonCalc = exports;
	} else {
		MoonCalc = global.MoonCalc = {};
	}


	// shortcuts for easier to read formulas

	var PI   = Math.PI,
	    rad  = PI / 180,
	    sin  = Math.sin,
	    cos  = Math.cos,
	    tan  = Math.tan,
	    asin = Math.asin,
	    atan = Math.atan2,
	    acos = Math.acos;


	// sun calculations are based on http://aa.quae.nl/en/reken/zonpositie.html formulas


	// date/time constants and conversions

	var dayMs = 1000 * 60 * 60 * 24,
	    J1970 = 2440588,
	    J2000 = 2451545;

	function toJulian(date) {
		return date.valueOf() / dayMs - 0.5 + J1970;
	}
	function fromJulian(j) {
		return new Date((j + 0.5 - J1970) * dayMs);
	}
	function toDays(date) {
		return toJulian(date) - J2000;
	}


	// general calculations for position

	var e = rad * 23.4397; // obliquity of the Earth

	function getRightAscension(l, b) {
		return atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l));
	}
	function getDeclination(l, b) {
		return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l));
	}
	function getAzimuth(H, phi, dec) {
		return atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi));
	}
	function getAltitude(H, phi, dec) {
		return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H));
	}
	function getSiderealTime(d, lw) {
		return rad * (280.16 + 360.9856235 * d) - lw;
	}


	// general sun calculations

	function getSolarMeanAnomaly(d) {
		return rad * (357.5291 + 0.98560028 * d);
	}
	function getEquationOfCenter(M) {
		return rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M));
	}
	function getEclipticLongitude(M, C) {
		var P = rad * 102.9372; // perihelion of the Earth
		return M + C + P + PI;
	}

	// sun times configuration (angle, morning name, evening name)

	var times = [
		[-0.83, 'sunrise',       'sunset'      ],
		[ -0.3, 'sunriseEnd',    'sunsetStart' ],
		[   -6, 'dawn',          'dusk'        ],
		[  -12, 'nauticalDawn',  'nauticalDusk'],
		[  -18, 'nightEnd',      'night'       ],
		[    6, 'goldenHourEnd', 'goldenHour'  ]
	];

	// adds a custom time to the times config

	MoonCalc.addTime = function (angle, riseName, setName) {
		times.push([angle, riseName, setName]);
	};


	// calculations for sun times

	var J0 = 0.0009;

	function getJulianCycle(d, lw) {
		return Math.round(d - J0 - lw / (2 * PI));
	}
	function getApproxTransit(Ht, lw, n) {
		return J0 + (Ht + lw) / (2 * PI) + n;
	}
	function getSolarTransitJ(ds, M, L) {
		return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L);
	}
	function getHourAngle(h, phi, d) {
		return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)));
	}

	// moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html formulas

	function getMoonCoords(d) { // geocentric ecliptic coordinates of the moon

		var L = rad * (218.316 + 13.176396 * d), // ecliptic longitude
		    M = rad * (134.963 + 13.064993 * d), // mean anomaly
		    F = rad * (93.272 + 13.229350 * d),  // mean distance

		    l  = L + rad * 6.289 * sin(M), // longitude
		    b  = rad * 5.128 * sin(F),     // latitude
		    dt = 385001 - 20905 * cos(M);  // distance to the moon in km

		return {
			ra: getRightAscension(l, b),
			dec: getDeclination(l, b),
			dist: dt
		};
	}

	MoonCalc.getMoonPosition = function (date, lat, lng) {

		var lw  = rad * -lng,
		    phi = rad * lat,
		    d   = toDays(date),

		    c = getMoonCoords(d),
		    H = getSiderealTime(d, lw) - c.ra,
		    h = getAltitude(H, phi, c.dec);

		// altitude correction for refraction
		h = h + rad * 0.017 / tan(h + rad * 10.26 / (h + rad * 5.10));

		return {
			azimuth: getAzimuth(H, phi, c.dec),
			altitude: h,
			distance: c.dist
		};
	};


	// calculations for illuminated fraction of the moon,
	// based on http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas

	MoonCalc.getMoonFraction = function (date) {

		var d = toDays(date),
		    s = getSunCoords(d),
		    m = getMoonCoords(d),

		    sdist = 149598000, // distance from Earth to Sun in km

		    phi = acos(sin(s.dec) * sin(m.dec) + cos(s.dec) * cos(m.dec) * cos(s.ra - m.ra)),
		    inc = atan(sdist * sin(phi), m.dist - sdist * cos(phi));

		return (1 + cos(inc)) / 2;
	};

	//////////////////////////////////////
	///////////////Get moon times/////////
	//////////////////////////////////////
	var Sky = [0.0, 0.0, 0.0];
	var RAn = [0.0, 0.0, 0.0];
	var Dec = [0.0, 0.0, 0.0];
	var VHz = [0.0, 0.0, 0.0];
	var Rise_time = [0, 0];
	var Set_time  = [0, 0];
	var Moonrise = false;
	var Moonset  = false;
	var DR = PI/180;
	var K1 = 15*DR*1.0027379;

	var Rise_az = 0.0;
	var Set_az  = 0.0;
	// determine Julian day from calendar date
	// (Jean Meeus, "Astronomical Algorithms", Willmann-Bell, 1991)
	function julian_day(date)
	{
	    var a, b, jd;
	    var gregorian;

	    var month = date.getMonth() + 1;
	    var day   = date.getDate();
	    var year  = date.getFullYear();

	    gregorian = (year < 1583) ? false : true;
	    
	    if ((month == 1)||(month == 2))
	    {
	        year  = year  - 1;
	        month = month + 12;
	    }

	    a = Math.floor(year/100);
	    if (gregorian) b = 2 - a + Math.floor(a/4);
	    else           b = 0.0;

	    jd = Math.floor(365.25*(year + 4716)) 
	       + Math.floor(30.6001*(month + 1)) 
	       + day + b - 1524.5;
	    
	    return jd;
	}

	// Local Sidereal Time for zone
	function lst( lon, jd, z )
	{
	    var s = 24110.5 + 8640184.812999999*jd/36525 + 86636.6*z + 86400*lon;
	    s = s/86400;
	    s = s - Math.floor(s);
	    return s*360*DR;
	}

	// moon's position using fundamental arguments 
	// (Van Flandern & Pulkkinen, 1979)
	function moon( jd )
	{
	    var d, f, g, h, m, n, s, u, v, w;

	    h = 0.606434 + 0.03660110129*jd;
	    m = 0.374897 + 0.03629164709*jd;
	    f = 0.259091 + 0.0367481952 *jd;
	    d = 0.827362 + 0.03386319198*jd;
	    n = 0.347343 - 0.00014709391*jd;
	    g = 0.993126 + 0.0027377785 *jd;

	    h = h - Math.floor(h);
	    m = m - Math.floor(m);
	    f = f - Math.floor(f);
	    d = d - Math.floor(d);
	    n = n - Math.floor(n);
	    g = g - Math.floor(g);

	    h = h*2*PI;
	    m = m*2*PI;
	    f = f*2*PI;
	    d = d*2*PI;
	    n = n*2*PI;
	    g = g*2*PI;

	    v = 0.39558*Math.sin(f + n);
	    v = v + 0.082  *Math.sin(f);
	    v = v + 0.03257*Math.sin(m - f - n);
	    v = v + 0.01092*Math.sin(m + f + n);
	    v = v + 0.00666*Math.sin(m - f);
	    v = v - 0.00644*Math.sin(m + f - 2*d + n);
	    v = v - 0.00331*Math.sin(f - 2*d + n);
	    v = v - 0.00304*Math.sin(f - 2*d);
	    v = v - 0.0024 *Math.sin(m - f - 2*d - n);
	    v = v + 0.00226*Math.sin(m + f);
	    v = v - 0.00108*Math.sin(m + f - 2*d);
	    v = v - 0.00079*Math.sin(f - n);
	    v = v + 0.00078*Math.sin(f + 2*d + n);
	    
	    u = 1 - 0.10828*Math.cos(m);
	    u = u - 0.0188 *Math.cos(m - 2*d);
	    u = u - 0.01479*Math.cos(2*d);
	    u = u + 0.00181*Math.cos(2*m - 2*d);
	    u = u - 0.00147*Math.cos(2*m);
	    u = u - 0.00105*Math.cos(2*d - g);
	    u = u - 0.00075*Math.cos(m - 2*d + g);
	    
	    w = 0.10478*Math.sin(m);
	    w = w - 0.04105*Math.sin(2*f + 2*n);
	    w = w - 0.0213 *Math.sin(m - 2*d);
	    w = w - 0.01779*Math.sin(2*f + n);
	    w = w + 0.01774*Math.sin(n);
	    w = w + 0.00987*Math.sin(2*d);
	    w = w - 0.00338*Math.sin(m - 2*f - 2*n);
	    w = w - 0.00309*Math.sin(g);
	    w = w - 0.0019 *Math.sin(2*f);
	    w = w - 0.00144*Math.sin(m + n);
	    w = w - 0.00144*Math.sin(m - 2*f - n);
	    w = w - 0.00113*Math.sin(m + 2*f + 2*n);
	    w = w - 0.00094*Math.sin(m - 2*d + g);
	    w = w - 0.00092*Math.sin(2*m - 2*d);

	    s = w/Math.sqrt(u - v*v);                  // compute moon's right ascension ...  
	    Sky[0] = h + Math.atan(s/Math.sqrt(1 - s*s));

	    s = v/Math.sqrt(u);                        // declination ...
	    Sky[1] = Math.atan(s/Math.sqrt(1 - s*s));

	    Sky[2] = 60.40974*Math.sqrt( u );          // and parallax
	}

	// 3-point interpolation
	function interpolate( f0, f1, f2, p )
	{
	    var a = f1 - f0;
	    var b = f2 - f1 - a;
	    var f = f0 + p*(2*a + b*(2*p - 1));

	    return f;
	}

	// returns value for sign of argument
	function sgn( x )
	{
	    var rv;
	    if (x > 0.0)      rv =  1;
	    else if (x < 0.0) rv = -1;
	    else              rv =  0;
	    return rv;
	}

	function test_moon( k, zone, t0, lat, plx )
	{
	    var ha = [0.0, 0.0, 0.0];
	    var a, b, c, d, e, s, z;
	    var hr, min, time;
	    var az, hz, nz, dz;

	    if (RAn[2] < RAn[0])
	        RAn[2] = RAn[2] + 2*PI;
	    
	    ha[0] = t0 - RAn[0] + k*K1;
	    ha[2] = t0 - RAn[2] + k*K1 + K1;
	    
	    ha[1]  = (ha[2] + ha[0])/2;                // hour angle at half hour
	    Dec[1] = (Dec[2] + Dec[0])/2;              // declination at half hour

	    s = Math.sin(DR*lat);
	    c = Math.cos(DR*lat);

	    // refraction + sun semidiameter at horizon + parallax correction
	    z = Math.cos(DR*(90.567 - 41.685/plx));

	    if (k <= 0)                                // first call of function
	        VHz[0] = s*Math.sin(Dec[0]) + c*Math.cos(Dec[0])*Math.cos(ha[0]) - z;

	    VHz[2] = s*Math.sin(Dec[2]) + c*Math.cos(Dec[2])*Math.cos(ha[2]) - z;
	    
	    if (sgn(VHz[0]) == sgn(VHz[2]))
	        return VHz[2];                         // no event this hour
	    
	    VHz[1] = s*Math.sin(Dec[1]) + c*Math.cos(Dec[1])*Math.cos(ha[1]) - z;

	    a = 2*VHz[2] - 4*VHz[1] + 2*VHz[0];
	    b = 4*VHz[1] - 3*VHz[0] - VHz[2];
	    d = b*b - 4*a*VHz[0];

	    if (d < 0)
	        return VHz[2];                         // no event this hour
	    
	    d = Math.sqrt(d);
	    e = (-b + d)/(2*a);

	    if (( e > 1 )||( e < 0 ))
	        e = (-b - d)/(2*a);

	    time = k + e + 1/120;                      // time of an event + round up
	    hr   = Math.floor(time);
	    min  = Math.floor((time - hr)*60);

	    hz = ha[0] + e*(ha[2] - ha[0]);            // azimuth of the moon at the event
	    nz = -Math.cos(Dec[1])*Math.sin(hz);
	    dz = c*Math.sin(Dec[1]) - s*Math.cos(Dec[1])*Math.cos(hz);
	    az = Math.atan2(nz, dz)/DR;
    	if (az < 0) az = az + 360;

	    if ((VHz[0] < 0)&&(VHz[2] > 0))
	    {
	        Rise_time[0] = hr;
	        Rise_time[1] = min;
	        Rise_az = az;
	        Moonrise = true;
	    }
	    
	    if ((VHz[0] > 0)&&(VHz[2] < 0))
	    {
	        Set_time[0] = hr;
	        Set_time[1] = min;
	        Set_az = az;
	        Moonset = true;
	    }

	    return VHz[2];
	}

	// format a positive integer with leading zeroes
	function zintstr( num, width )
	{
	    var str = num.toString(10);
	    var len = str.length;
	    var intgr = "";
	    var i;

	    for (i = 0; i < width - len; i++)          // append leading zeroes
	        intgr += '0';

	    for (i = 0; i < len; i++)                  // append digits
	        intgr += str.charAt(i);

	    return intgr;
	}

	// format a real number
	function frealstr( num, width, fract )
	{
	    var str = num.toFixed(fract);
	    var len = str.length;
	    var real = "";
	    var i;

	    for (i = 0; i < width - len; i++)          // append leading spaces
	        real += ' ';

	    for (i = 0; i < len; i++)                  // append digits
	        real += str.charAt(i);

	    return real;
	}

	MoonCalc.getMoonTimes = function riseset( lat, lon ,date)
	{
	    var i, j, k;
	    var zone = Math.round(date.getTimezoneOffset()/60);
	    var jd = julian_day(date) - 2451545;           // Julian day relative to Jan 1.5, 2000
	    
	//    if ((sgn(zone) == sgn(lon))&&(zone != 0))
	//        window.alert("WARNING: time zone and longitude are incompatible!");

	    var mp = new Array(3);                     // create a 3x3 array
	    for (i = 0; i < 3; i++)
	    {
	        mp[i] = new Array(3);
	        for (j = 0; j < 3; j++)
	            mp[i][j] = 0.0;
	    }

	/////////// new stuff

	    var x = lon;
	    zone = Math.round(-x/15);

	///////////

	    lon = lon/360;
	    var tz = zone/24;
	    var t0 = lst(lon, jd, tz);                 // local sidereal time

	    jd = jd + tz;                              // get moon position at start of day

	    for (k = 0; k < 3; k++)
	    {
	        moon(jd);
	        mp[k][0] = Sky[0];
	        mp[k][1] = Sky[1];
	        mp[k][2] = Sky[2];
	        jd = jd + 0.5;      
	    }   

	    if (mp[1][0] <= mp[0][0])
	        mp[1][0] = mp[1][0] + 2*PI;

	    if (mp[2][0] <= mp[1][0])
	        mp[2][0] = mp[2][0] + 2*PI;

	    RAn[0] = mp[0][0];
	    Dec[0] = mp[0][1];

	    Moonrise = false;                          // initialize
	    Moonset  = false;
	    
	    for (k = 0; k < 24; k++)                   // check each hour of this day
	    {
	        var ph = (k + 1)/24;
	        
	        RAn[2] = interpolate(mp[0][0], mp[1][0], mp[2][0], ph);
	        Dec[2] = interpolate(mp[0][1], mp[1][1], mp[2][1], ph);
	        
	        VHz[2] = test_moon(k, zone, t0, lat, mp[1][2]);

	        RAn[0] = RAn[2];                       // advance to next hour
	        Dec[0] = Dec[2];
	        VHz[0] = VHz[2];
	    }

	    // return results
	    var moonsetTime;
	    var moonriseTime;
	    var nextDay =false;
	    if (zintstr( Set_time[0], 2)<zintstr(Rise_time[0], 2))
	    {
	    	moonsetTime = new Date(date.getFullYear(), date.getMonth(), date.getDate()+1,zintstr( Set_time[0], 2), zintstr( Set_time[1], 2));
	    	moonriseTime = new Date(date.getFullYear(), date.getMonth(), date.getDate(), zintstr(Rise_time[0], 2), zintstr(Rise_time[1], 2));
	    	nextDay=true;
	    }
	    else{
	    	moonsetTime = new Date(date.getFullYear(), date.getMonth(), date.getDate(),zintstr( Set_time[0], 2), zintstr( Set_time[1], 2));
	    	moonriseTime = new Date(date.getFullYear(), date.getMonth(), date.getDate(), zintstr(Rise_time[0], 2), zintstr(Rise_time[1], 2));
	    	nextDay=false;
	    }
	    
	    var moonTimes=[moonriseTime,moonsetTime,nextDay];
	    return moonTimes;
	}

}(this));