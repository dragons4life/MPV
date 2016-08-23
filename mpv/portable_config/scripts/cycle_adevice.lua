--binds a hotkey to cycle through available audio devices
local numdevices = 0
local current = false
local currentd = ""
local devicenames = {};
local devicedesc = {};
local devices = {}
local goBack = false

local function cycle_adevice()
	numdevices = 0
	deviceList = mp.get_property_native("audio-device-list")

	for index, e in ipairs(deviceList) do
		if string.find(e.name, "wasapi", 1, true) then
			numdevices = numdevices + 1
			devicenames[numdevices] = e.name;
			devicedesc[numdevices] = e.description;
		end

    end

	for i=1, numdevices do
		currentd = devicenames[i]

		if string.find(mp.get_property("audio-device"), currentd, 1, true) or string.find(mp.get_property("audio-device"), "auto", 1, true) then
			current = true
			if goBack then
				i = i - 2
				goBack = false
				if i == -1 then i = numdevices - 1; end
			end

			if i == numdevices then i = 0; end
		end

		if current then
			mp.set_property("audio-device",devicenames[i+1])
			print("audio="..devicedesc[i+1])
			mp.osd_message("audio="..devicedesc[i+1])
			current = false
			break
		end

	end

end

local function cycle_back()
	goBack = true
	cycle_adevice()
end


mp.add_key_binding("a", "cycle_adevice", cycle_adevice)
mp.add_key_binding("A", "cycleBack_adevice", cycle_back)
