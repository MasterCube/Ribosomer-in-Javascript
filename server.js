// server.js
console.log("server is starting...");
// import the package express:
var express = require('express');
// this all package is actually a function that
// I stored in the express variable:
var app = express();
var server = app.listen(3000, listening);
function listening(){
  console.log('I am listening...');
}

// NEW - Add CORS headers - see https://enable-cors.org/server_expressjs.html
app.use(function(req, res, next) {
  res.header("Access-Control-Allow-Origin", "*");
  res.header(
    "Access-Control-Allow-Headers",
    "Origin, X-Requested-With, Content-Type, Accept"
  );
  next();
});
// Now I can host static files
// app.use(express.static('public'));
// use express capabilities to host html
// files, jpeg files, movie files, ...
//app.use(express.static('website'));
app.use(express.static('j5sketch RIBOSOMER'));

// let allowCrossDomain = function(req, res, next) {
//   res.header('Access-Control-Allow-Origin', "*");
//   res.header('Access-Control-Allow-Headers', "*");
//   next();
// }
// app.use(allowCrossDomain);

// setting up a route. The purpose is to send requests to the server
// and get response from it. The request is written as a route (a path)
// to the 'public=website' directory which is 'localhost:3000'. In the
// public directory you have created a route (path) like /flower.
// Upon calling flower, a particular response is designed. See example:
app.get('/flower', sendFlower);
// // the first arg is the request,
// // the second is the callback and will be the response:
  function sendFlower(request, response){
  response.send("my response is that I love flowers too. Jouarez is the best!");
}

// next step:
app.get('/search/:flower', sendFlower);
// the first arg is the request,
// the second is the callback and will be the response:
function sendFlower(request, response){
  var data = request.params;
  response.send("my response is now that I love " + data.flower +  " too. Jouarez is still the best!");
}
