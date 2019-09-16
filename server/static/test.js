let window = {
  a: 1,
  b: () => {}
}

let w = Object.assign({}, window)

for (let i in window) {
  try {
      Object.defineProperty(window, i, {
          get: function () {
              console.log('name: ', i)
              console.log('value: ', w[i])
              return w[i]
          }
      })
  } catch(err){}
}

console.log(window.a)
console.log(window.b)
