/* eslint-disable @typescript-eslint/no-explicit-any */
export class QWebChannel {
    transport: any;
    execCallbacks:{[key: number]: (data: any) => void};
    execId: number;
    objects: {[key: string]: any};

    constructor(transport: any, initCallback: (channel: QWebChannel) => void) {
        this.transport = transport;
        this.execCallbacks = {};
        this.execId = 0;
        this.objects = {};

        this.transport.onmessage = (message: any) => {
            const data = JSON.parse(message.data);
            if (data.type == 1){                //prop update
                for (const  i in data.properties){
                    const prop = data.properties[i];
                    this.objects[prop.object].wrappedProperties[prop.property] = prop.value;
                }
            } else if (data.type == 2){         //response
                const callback = this.execCallbacks[data.id];
                if (callback){
                    delete this.execCallbacks[data.id];
                    callback(data.response);
                }
            }
        };
        this.transport.send(JSON.stringify({type: 0})); //Init connection   
        initCallback(this);                             //Expose backend hook
    }
}