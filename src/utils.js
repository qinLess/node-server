const readline = require('readline')
const fs = require('fs')
const superagent = require('superagent')

const cats = '[{"count":3561,"id":-100,"ids":[-100,1,7,8,2,4,5,6,3],"level":1,"name":"美食","sub_categories":[{"count":3561,"id":-100,"image_url":"8644e1323e3a11137d7478d2c83095d4png","level":1,"name":"全部"},{"count":970,"id":1,"image_url":"d8515995a32586701f3c184f422cf66apng","level":1,"name":"简餐便当"},{"count":850,"id":7,"image_url":"18eedac58bbcfd59e4a2805abfe11ff2jpeg","level":1,"name":"小吃炸串"},{"count":722,"id":8,"image_url":"6f8c0289b217b04fffd70c5f1f397103png","level":1,"name":"地方菜系"},{"count":656,"id":2,"image_url":"af571d2ac336d2bde99a2f2b7be97505png","level":1,"name":"面食粥点"},{"count":136,"id":4,"image_url":"f36695c9aab9774acbf5ba7e553b8a99jpeg","level":1,"name":"香锅冒菜"},{"count":91,"id":5,"image_url":"ce8b3161428d69729cdeb8095946326bjpeg","level":1,"name":"日韩料理"},{"count":78,"id":6,"image_url":"a5766bbe0b3ce5c74d81b63e06b614e1png","level":1,"name":"轻食西餐"},{"count":58,"id":3,"image_url":"808fc83b1b9fe47e04219466b269c3bapng","level":1,"name":"汉堡披萨"}]},{"count":1296,"id":207,"ids":[207],"level":1,"name":"快餐便当","sub_categories":[{"count":1296,"id":207,"image_url":"44545a0518aab93817cfe611e88bb702png","level":1,"name":"全部"},{"count":599,"id":265,"image_url":"be84bc4d7cf12deee9115b16eb099302png","level":2,"name":"简餐"},{"count":379,"id":213,"image_url":"02e6c9e3bf338ec0ba0d923717b9f8acpng","level":2,"name":"米粉面馆"},{"count":293,"id":209,"image_url":"66b78c0e7099c278977298d7c6042c80png","level":2,"name":"盖浇饭"},{"count":131,"id":215,"image_url":"af6ab89041b3e77fe115d1e4b72d69f0png","level":2,"name":"包子粥店"},{"count":106,"id":217,"image_url":"65b575c2278a3f6e5c70af45b578cbeepng","level":2,"name":"饺子馄饨"},{"count":82,"id":214,"image_url":"4d347d0dc65dd75fb2911256aabf2327png","level":2,"name":"麻辣烫"},{"count":54,"id":219,"image_url":"eddd9dc7e5d21debe2fb278ae01fefe9png","level":2,"name":"香锅砂锅"},{"count":40,"id":216,"image_url":"4c6af48f68284ad91c6d95d2bd3f4aa6png","level":2,"name":"生煎锅贴"},{"count":37,"id":267,"image_url":"79637dc36d67de4fe48d121ea77b3eddpng","level":2,"name":"黄焖鸡米饭"},{"count":29,"id":212,"image_url":"7d47af01fccc46fc3621865a9cc07c93png","level":2,"name":"汉堡"},{"count":27,"id":269,"image_url":"180cb951c2d4eb2e220debf4571bf83apng","level":2,"name":"煲仔饭"},{"count":8,"id":266,"image_url":"c09d1ff71384e2e1664f72e0a928810dpng","level":2,"name":"烧腊饭"},{"count":6,"id":268,"image_url":"cdf208b399b854e456f23d28b1972e97png","level":2,"name":"咖喱饭"}]},{"count":611,"id":220,"ids":[220],"level":1,"name":"特色菜系","sub_categories":[{"count":611,"id":220,"image_url":"6f8c0289b217b04fffd70c5f1f397103png","level":1,"name":"全部"},{"count":236,"id":225,"image_url":"2d098842683548f9626cf0a8c879257dpng","level":2,"name":"江浙菜"},{"count":149,"id":263,"image_url":"94ac841e2c3e27f8eeeaa917574ed574png","level":2,"name":"其他菜系"},{"count":140,"id":221,"image_url":"43b0e4694f8ebc393cce6723d5df5222png","level":2,"name":"川湘菜"},{"count":99,"id":232,"image_url":"a33f1ec0044ddd4d282fbc8b1f0a946fpng","level":2,"name":"海鲜"},{"count":52,"id":231,"image_url":"c03d81f550eb849ed2d4d0290ced9099png","level":2,"name":"火锅烤鱼"},{"count":27,"id":222,"image_url":"e320bf1ab9762cb1faad27d79f51219cpng","level":2,"name":"粤菜"},{"count":8,"id":223,"image_url":"aa4de1e9b54170cf495d8052407658c5png","level":2,"name":"东北菜"},{"count":5,"id":226,"image_url":"741d15270496d7699dd2e7804fccc7a1png","level":2,"name":"西北菜"},{"count":3,"id":228,"image_url":"a7e6d9cf1993fa4fe0bd02d74d40c9c2png","level":2,"name":"新疆菜"},{"count":2,"id":224,"image_url":"54dabf93116f4a336fcc91431be43828png","level":2,"name":"云南菜"},{"count":1,"id":227,"image_url":"e19bf59188a157dfc372b3d254fc986dpng","level":2,"name":"鲁菜"}]},{"count":176,"id":260,"ids":[260],"level":1,"name":"异国料理","sub_categories":[{"count":176,"id":260,"image_url":"754c5c2ad1b01668a7186ec5f0fb0e59png","level":1,"name":"全部"},{"count":91,"id":229,"image_url":"cf8c84a2fe5ecf27b21bcbddc1724d36png","level":2,"name":"日韩料理"},{"count":69,"id":230,"image_url":"78c45200d58e5c02cb70fb8287df732dpng","level":2,"name":"西餐"},{"count":29,"id":211,"image_url":"bb7eb2afe778ba9afbe54f9d282818d1png","level":2,"name":"披萨意面"},{"count":9,"id":264,"image_url":"614053401fddc171eed0436f3cd1f7dcpng","level":2,"name":"东南亚菜"}]},{"count":850,"id":233,"ids":[233],"level":1,"name":"小吃夜宵","sub_categories":[{"count":850,"id":233,"image_url":"7d714540b1590552d991fd731e8772a3png","level":1,"name":"全部"},{"count":410,"id":237,"image_url":"90483b16d9598aec798263220eb3a821png","level":2,"name":"地方小吃"},{"count":245,"id":236,"image_url":"d049fb617edcea921185258d1675db83png","level":2,"name":"小龙虾"},{"count":125,"id":234,"image_url":"71164ef684e8a13b5e66a20a1c55671cpng","level":2,"name":"炸鸡炸串"},{"count":115,"id":218,"image_url":"3c6e2763cf4ee56f18fd1b7360585fb3png","level":2,"name":"烧烤"},{"count":68,"id":235,"image_url":"efdba78945f83ed1e8e6e838718b4c65png","level":2,"name":"鸭脖卤味"},{"count":30,"id":238,"image_url":"d7e0be7e5420e213ea42e4fa3efa762bpng","level":2,"name":"零食"}]},{"count":911,"id":-102,"ids":[-102,11,9,12,10],"level":1,"name":"甜品饮品","sub_categories":[{"count":911,"id":-102,"image_url":"cac2f06034f8e73eb7793b08c1987049png","level":1,"name":"全部"},{"count":364,"id":11,"image_url":"62c31de0a7f41231e4c6934dea621a33png","level":1,"name":"奶茶果汁"},{"count":294,"id":9,"image_url":"213cbac0242d4845d1d28af0fa5fe35epng","level":1,"name":"甜品"},{"count":183,"id":12,"image_url":"ac94b005c97ef158282326cb49389893png","level":1,"name":"面包蛋糕"},{"count":70,"id":10,"image_url":"c2f05ef82a7ee44b7848b7fb598d42e3png","level":1,"name":"咖啡"}]},{"count":233,"id":244,"ids":[244],"level":1,"name":"果蔬生鲜","sub_categories":[{"count":233,"id":244,"image_url":"1ce198f37a81285f4afa2aaf826a558fpng","level":1,"name":"全部"},{"count":157,"id":245,"image_url":"a831a37ec670ca93cd35a8a6b5a20e62png","level":2,"name":"水果"},{"count":65,"id":246,"image_url":"1729548b88614c1b3a6e71ef7f89f294png","level":2,"name":"蔬菜豆品"},{"count":61,"id":247,"image_url":"6d3cef77e055d03598cba821ebcf1f06png","level":2,"name":"肉禽蛋品"},{"count":13,"id":270,"image_url":"a2ab438ee4ac09e6e53b3f96694bac81png","level":2,"name":"海鲜水产"}]},{"count":181,"id":252,"ids":[252],"level":1,"name":"商店超市","sub_categories":[{"count":181,"id":252,"image_url":"df21b511f287ccb402e68285d2653caepng","level":1,"name":"全部"},{"count":111,"id":254,"image_url":"92ae70438be9a3adfc5a560c1e6ae818png","level":2,"name":"大型超市"},{"count":48,"id":271,"image_url":"841d136b17fa4cb871a296c9e4997cfapng","level":2,"name":"便利店"},{"count":42,"id":273,"image_url":"c2b0e2b27ea55a9a7211f14ad95dcd0apng","level":2,"name":"休闲零食"},{"count":9,"id":274,"image_url":"7df84232aebbb5ffb53e564c9e328d31png","level":2,"name":"名酒坊"},{"count":2,"id":255,"image_url":"825031dc99e1f99c26feb7186b6cf3a6png","level":2,"name":"水站"},{"count":1,"id":257,"image_url":"b435af6662fd0b3e9fb6537474753f72png","level":2,"name":"粮油副食"}]},{"count":249,"id":275,"ids":[275],"level":1,"name":"鲜花绿植","sub_categories":[{"count":249,"id":251,"image_url":"cf598de7338b4bf9dd2924736c4ec9d2png","level":2,"name":"浪漫鲜花"}]},{"count":13,"id":276,"ids":[276],"level":1,"name":"医药健康","sub_categories":[{"count":13,"id":277,"image_url":"616af0303ec6775365631b7aa1df6da0png","level":2,"name":"医药健康"}]},{"count":1075,"id":-103,"ids":[-103,-5,-9,-6,-8,-7,-4,-11],"level":1,"name":"早餐","sub_categories":[{"count":1075,"id":-103,"image_url":"f2c5870ff4c3f9fd702346ae907ab038png","level":1,"name":"全部"},{"count":379,"id":-5,"image_url":"02e6c9e3bf338ec0ba0d923717b9f8acpng","level":1,"name":"米粉面馆"},{"count":364,"id":-9,"image_url":"62c31de0a7f41231e4c6934dea621a33png","level":1,"name":"奶茶果汁"},{"count":131,"id":-6,"image_url":"af6ab89041b3e77fe115d1e4b72d69f0png","level":1,"name":"包子粥店"},{"count":106,"id":-8,"image_url":"65b575c2278a3f6e5c70af45b578cbeepng","level":1,"name":"饺子馄饨"},{"count":40,"id":-7,"image_url":"4c6af48f68284ad91c6d95d2bd3f4aa6png","level":1,"name":"生煎锅贴"},{"count":29,"id":-4,"image_url":"7d47af01fccc46fc3621865a9cc07c93png","level":1,"name":"汉堡"},{"count":26,"id":-11,"image_url":"512232422a83e25a2c0a5588b7b6e730png","level":1,"name":"面包"}]},{"count":1052,"id":-104,"ids":[-104,-15,-12,-14,-13],"level":1,"name":"午餐","sub_categories":[{"count":1052,"id":-104,"image_url":"b0c0d7e523251cec33b96953d72c0348jpeg","level":1,"name":"全部"},{"count":599,"id":-15,"image_url":"be84bc4d7cf12deee9115b16eb099302png","level":1,"name":"简餐"},{"count":293,"id":-12,"image_url":"66b78c0e7099c278977298d7c6042c80png","level":1,"name":"盖浇饭"},{"count":131,"id":-14,"image_url":"af6ab89041b3e77fe115d1e4b72d69f0png","level":1,"name":"包子粥店"},{"count":29,"id":-13,"image_url":"7d47af01fccc46fc3621865a9cc07c93png","level":1,"name":"汉堡"}]},{"count":1671,"id":-105,"ids":[-105,-21,-16,-17,-20,-19,-18,-22],"level":1,"name":"下午茶","sub_categories":[{"count":1671,"id":-105,"image_url":"eb63a10e117d6ec35d6ceed2df72eadapng","level":1,"name":"全部"},{"count":535,"id":-21,"image_url":"71164ef684e8a13b5e66a20a1c55671cpng","level":1,"name":"炸鸡小吃"},{"count":364,"id":-16,"image_url":"62c31de0a7f41231e4c6934dea621a33png","level":1,"name":"奶茶果汁"},{"count":294,"id":-17,"image_url":"213cbac0242d4845d1d28af0fa5fe35epng","level":1,"name":"甜品"},{"count":183,"id":-20,"image_url":"ac94b005c97ef158282326cb49389893png","level":1,"name":"面包蛋糕"},{"count":157,"id":-19,"image_url":"a831a37ec670ca93cd35a8a6b5a20e62png","level":1,"name":"水果"},{"count":70,"id":-18,"image_url":"c2f05ef82a7ee44b7848b7fb598d42e3png","level":1,"name":"咖啡"},{"count":68,"id":-22,"image_url":"efdba78945f83ed1e8e6e838718b4c65png","level":1,"name":"鸭脖卤味"}]},{"count":2450,"id":-107,"ids":[-107,-31,-35,-33,-37,-32,-34,-38,-36,-40,-39],"level":1,"name":"晚餐","sub_categories":[{"count":2450,"id":-107,"image_url":"07af2e54dddf2d82a05ecd1f2a194e13png","level":1,"name":"全部"},{"count":599,"id":-31,"image_url":"be84bc4d7cf12deee9115b16eb099302png","level":1,"name":"简餐"},{"count":410,"id":-35,"image_url":"90483b16d9598aec798263220eb3a821png","level":1,"name":"地方小吃"},{"count":379,"id":-33,"image_url":"02e6c9e3bf338ec0ba0d923717b9f8acpng","level":1,"name":"米粉面馆"},{"count":364,"id":-37,"image_url":"62c31de0a7f41231e4c6934dea621a33png","level":1,"name":"奶茶果汁"},{"count":293,"id":-32,"image_url":"66b78c0e7099c278977298d7c6042c80png","level":1,"name":"盖浇饭"},{"count":140,"id":-34,"image_url":"43b0e4694f8ebc393cce6723d5df5222png","level":1,"name":"川湘菜"},{"count":125,"id":-38,"image_url":"71164ef684e8a13b5e66a20a1c55671cpng","level":1,"name":"炸鸡炸串"},{"count":82,"id":-36,"image_url":"4d347d0dc65dd75fb2911256aabf2327png","level":1,"name":"麻辣烫"},{"count":29,"id":-40,"image_url":"bb7eb2afe778ba9afbe54f9d282818d1png","level":1,"name":"披萨意面"},{"count":29,"id":-39,"image_url":"7d47af01fccc46fc3621865a9cc07c93png","level":1,"name":"汉堡"}]},{"count":1385,"id":-106,"ids":[-106,-29,-24,-28,-27,-26,-25,-23],"level":1,"name":"夜宵","sub_categories":[{"count":1385,"id":-106,"image_url":"d763d496b0454edd5ddb7eb5d796b56bpng","level":1,"name":"全部"},{"count":410,"id":-29,"image_url":"90483b16d9598aec798263220eb3a821png","level":1,"name":"地方小吃"},{"count":379,"id":-24,"image_url":"02e6c9e3bf338ec0ba0d923717b9f8acpng","level":1,"name":"米粉面馆"},{"count":245,"id":-28,"image_url":"d049fb617edcea921185258d1675db83png","level":1,"name":"小龙虾"},{"count":125,"id":-27,"image_url":"71164ef684e8a13b5e66a20a1c55671cpng","level":1,"name":"炸鸡炸串"},{"count":115,"id":-26,"image_url":"3c6e2763cf4ee56f18fd1b7360585fb3png","level":1,"name":"烧烤"},{"count":82,"id":-25,"image_url":"4d347d0dc65dd75fb2911256aabf2327png","level":1,"name":"麻辣烫"},{"count":29,"id":-23,"image_url":"7d47af01fccc46fc3621865a9cc07c93png","level":1,"name":"汉堡"}]}]'
const ak = 'E4805d16520de693a3fe707cdc962045'

function getParentCat(cat) {
  var pcats = [];
  var json = JSON.parse(cats);
  for (var i in json) {
    var sub_categories = json[i].sub_categories;
    if (sub_categories) {
      for (var j in sub_categories) {
        if (sub_categories[j].name == cat) {
          pcats.push(json[i].name);
          break
        }
      }
    }
  }
  if (pcats.length == 0) {
    pcats.push("其他")
  }
  return pcats
}

async function getLocation(lat, lng) {
  let url = `http://api.map.baidu.com/geocoder?location=${lat},${lng}&output=json&key=${ak}`
  let ua = 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/72.0.3626.96 Safari/537.36'
  let res = await superagent.get(url).set('User-Agent', ua)
  try {
    return JSON.parse(res.text)
  } catch (err) {
    console.log(res)
    process.exit()
  }
}

function getImagePath(t) {
  if (t) {
    var e = {
      supportedWebp: false
    }
    var n = "//fuss10.elemecdn.com"
    var a = t.match(/\|(\d+)x(\d+)/)
    var r = /\/dist|static\.alpha\.elenet\.me|static\.beta\.elenet\.me|static\.elenet\.me|static11\.elemecdn\.com|faas\.elemecdn\.com|cdn\.faas\.elenet\.me/.test(t);
    if (r) return t.replace(/\|(\d+)x(\d+)/, "");
    if (t = n + t.replace("http:", "").replace(n, "").replace(/^(\w)(\w\w)(\w{29}(\w*))(.+)?$/, "/$1/$2/$3.$4"), !a) return "" + t + (e.supportedWebp ? "?imageMogr2/format/webp/quality/85" : "")
    var i = t.split("|")[0],
      o = a[1],
      s = a[2],
      c = e.isRetina ? 2 : 1;
    return i += "?imageMogr2/thumbnail/" + o * c + "x" + s * c, e.supportedWebp && (i += "/format/webp/quality/85"), i
  }
  return "";
};


exports.readSyncByRl = async (tips) => {
  tips = tips || '> ';
  const rl = readline.createInterface({
    input: process.stdin,
    output: process.stdout
  })
  let answer = await rl.question(tips)
  rl.close()
  return answer.trim()
}

exports.save = (path, data) => {
  if (fs.existsSync(path)) {
    fs.unlinkSync(path);
  }
  fs.writeFileSync(path, data, 'utf-8');
}

// 解析数据
exports.paserData = async (dataList, dataObj) => {
  dataObj = dataObj || {}
  for (let dataItem of dataList) {
    dataItem = dataItem.restaurant
    let data = dataObj[dataItem.id]
    if (!data) {
      let categories = []
      let mainCate = []
      let flavors = dataItem.flavors || []
      for (let flavor of flavors) {
        let newFlavor = flavor
        let parentNames = getParentCat(newFlavor.name)
        mainCate = mainCate.concat(parentNames)
        categories.push(newFlavor.name)
      }
      let location = await getLocation(dataItem.latitude, dataItem.longitude)
      dataObj[dataItem.id] = {
        url: `https://www.ele.me/shop/${dataItem.id}`,
        shopId: dataItem.id,
        name: dataItem.name,
        categories: categories,
        mainCate: Array.from(new Set(mainCate)),

        province: location.result.addressComponent.province,
        city: location.result.addressComponent.city,

        logo: getImagePath(dataItem.image_path),
        latitude: dataItem.latitude,
        longitude: dataItem.longitude,
        score: dataItem.rating,
        isNew: dataItem.is_new > -1 ? '是' : '否',
        isPremium: dataItem.is_premium > -1 ? '是' : '否',
        phone: dataItem.phone,
        openingHours: dataItem.opening_hours,
        monthSalesCount: dataItem.recent_order_num,
        description: dataItem.piecewise_agent_fee.description,
        tips: dataItem.piecewise_agent_fee.tips,
        promotions: dataItem.activities,
        delivery: dataItem.delivery_mode ? dataItem.delivery_mode.text : '',
        avgDeliveryDuration: dataItem.order_lead_time,
        deliveryStartPrice: dataItem.float_minimum_order_amount,
        deliveryFee: dataItem.float_delivery_fee,
        supports: dataItem.supports
      }
    }
  }
  return dataObj
}

// 解析店铺商品
exports.parseProd = async (dataList) => {
  let productList = []
  for (let dataItem in dataList) {
    let foods = []
    for (let foodItem in dataItem.foods) {
      let skuDetail = []
      let price = null
      for (let skuItem in foodItem.specfoods) {
        skuDetail.push({
          sku: skuItem.name,
          price: skuItem.price,
          original_price: skuItem.original_price,
          stock: skuItem.stock,
          packingFee: skuItem.packing_fee
        })
        if (!price) {
          price = skuItem.price || skuItem.original_price  
        }
      }
      let imageUrl = foodItem.image_path ? image_base_url + foodItem.image_path.substring(0, 1) + '/' + foodItem.image_path.substring(1, 3) + '/' + foodItem.image_path.substring(3) + '.jpeg' : null

      foods.push({
        name: foodItem.name,
        monthSalesCount: foodItem.month_sales,
        score: foodItem.rating,
        pid: foodItem.virtual_food_id,
        desc: foodItem.description,
        price: price,
        images: imageUrl,
        recommendsCount: foodItem.satisfy_count,
        recommendsRate: foodItem.satisfy_rate,
        commentsCount: foodItem.rating_count,
        skuDetail: skuDetail
      })
    }
    productList.push({
      name: dataItem.name,
      products: foods
    })
  }
  return productList
}